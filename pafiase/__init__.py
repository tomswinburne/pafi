"""
pafiase — PAFI (Projected Average Force Integrator) with ASE calculators.

Implements the PAFI method [Swinburne & Marinica, PRL 2018] for computing
free-energy differences along a reference path using constrained Langevin MD.

A cubic-spline pathway is built from a list of NEB images. At each reaction
coordinate value ``r ∈ [0, 1]`` the system is constrained to move
perpendicular to the path tangent, and the projected force ``dF/dr`` is
accumulated. Integrating ``<dF/dr>`` over ``r`` gives the free-energy
profile ``F(r, T)``.

References:
    T. D. Swinburne and S. L. Marinica,
    "Unsupervised Calculation of Free Energy Barriers in Large Crystalline
    Systems", PRL 120, 135503 (2018).

(c) 2025 Thomas Swinburne — tomswinburne.github.io
"""

import json
import time
import uuid


import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline

from ase import units
from ase.constraints import FixedMode
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import FIRE2

def _naive_block_variance(all_data,barcodes=None,block_size=500):
    """
    Naive blocking analysis, assuming large minimal block
    pyblock tends to underestimate block size
    """
    if barcodes is None:
        n_blocks = max(1,all_data.size // block_size)
        return all_data.var() / float(n_blocks)
    else:
        _var = 0.0
        n_blocks = 0
        bcs = np.unique(barcodes)
        for bc in bcs:
            data = all_data[barcodes==bc]
            _var += data.var() / len(bcs)
            n_blocks += max(1,data.size // block_size)
        return _var / float(n_blocks)




    for bc in np.unique(barcodes):
            _vars.append(_single_block_variance(all_data[barcodes==bc]))
        
    


class PAFIData:
    """Container for PAFI simulation data with I/O and plotting.

    Can be used standalone to load and analyse results without a calculator::

        data = PAFIData()
        data.read_json("pafi_results.json")
        data.plot_data(filename="free_energy.png")
    """

    def __init__(self):
        self.reset_data()

    def reset_data(self):
        """Clear all accumulated simulation data."""
        self.data = {
            "barcode": [],
            "step": [],
            "r": [],
            "T": [],
            "target_T": [],
            "strain": [],
            "dF": [],
            "Psi": [],
            "dx.t": [],
            "max(dx)": [],
            "dt": [],
        }

    def write_json(self, filename):
        """Serialise ``self.data`` to a JSON file.

        Args:
            filename: Path to the output file (created or overwritten).
        """
        with open(filename, "w") as f:
            json.dump(self.data, f, indent=4)
        print(f"Written to {filename}")

    def read_json(self, filename):
        """Load simulation data from a JSON file into ``self.data``.

        Replaces any previously accumulated data.

        Args:
            filename: Path to a JSON file previously written by
                :meth:`write_json`.
        """
        with open(filename, "r") as f:
            self.data = json.load(f)
        print(f"Read from {filename}")

    def _strain_mask(self, data, strain):
        """Return a boolean mask selecting samples matching ``strain``.

        Args:
            data: Data dictionary.
            strain: Voigt 6-tuple/list to match, or ``None`` to select all.
        """
        if strain is None or "strain" not in data or len(data["strain"]) == 0:
            return np.ones(len(data["r"]), dtype=bool)
        strain = tuple(np.round(strain, 10))
        return np.array([tuple(np.round(s, 10)) == strain for s in data["strain"]])

    def _unique_strains(self, data):
        """Return sorted list of unique strain Voigt tuples in ``data``."""
        if "strain" not in data or len(data["strain"]) == 0:
            return [None]
        seen = set()
        result = []
        for s in data["strain"]:
            key = tuple(np.round(s, 10))
            if key not in seen:
                seen.add(key)
                result.append(key)
        return sorted(result)

    def compute_dF_statistics(self, data, T, strain=None):
        """Compute mean and variance of the projected force at each reaction coordinate.

        Selects samples from ``data`` matching temperature ``T`` and
        optionally ``strain``, then for each unique ``r`` value computes
        the block-averaged mean and variance of ``dF``.

        Args:
            data: Data dictionary in the format of ``self.data``.
            T: Target temperature (K).
            strain: Voigt 6-tuple to filter by, or ``None`` for all strains.

        Returns:
            Tuple ``(unique_r, dF_ave, dF_var)``.
        """
        int_T = np.array(data["target_T"]).astype(int)
        mask = (int_T == int(T)) & self._strain_mask(data, strain)
        sub = {k: np.array(data[k])[mask] for k in data}

        unique_r = np.sort(np.unique(sub["r"]))
        assert len(unique_r) > 1, (
            f"Only one unique value of 'r' found for T={T}; cannot compute statistics."
        )

        dF_ave, dF_var = [], []
        for r_val in unique_r:
            dF_data = sub["dF"][sub["r"] == r_val]
            bc_data = sub["barcode"][sub["r"] == r_val]
            dF_ave.append(-np.mean(dF_data))
            dF_var.append(_naive_block_variance(dF_data, bc_data))

        return unique_r, np.array(dF_ave), np.array(dF_var)

    def plot_data(self, data=None, filename=None, symmetric=False):
        """Plot the free-energy profile ``F(r)`` from accumulated PAFI data.

        Groups data by temperature, computes mean and variance of ``dF`` at
        each unique ``r`` via :meth:`compute_dF_statistics`, fits a spline,
        and integrates to obtain ``F(r)``.  Uncertainty bands (±1 σ) are
        shown as shaded regions.

        Args:
            data: Data dictionary in the format of ``self.data``.  If
                ``None``, uses ``self.data``.
            filename: If given, save the figure to this path; otherwise
                call ``plt.show()``.
            symmetric: Averages over forward and backward integration.
                Default False

        Returns:
            Tuple ``(r_dense, F_dense)`` — the dense reaction-coordinate
            array and the integrated free-energy profile for the last
            temperature group plotted, or ``None`` if plotting was skipped.

        Raises:
            AssertionError: If required keys ``r``, ``T``, or ``dF`` are
                absent from ``data``.
        """
        if data is None:
            data = self.data
        assert "r" in data, "Data must contain 'r' key"
        assert "T" in data, "Data must contain 'T' key"
        assert "dF" in data, "Data must contain 'dF' key"

        int_T = np.array(data["target_T"]).astype(int)
        unique_T = np.sort(np.unique(int_T))
        unique_strains = self._unique_strains(data)

        plt.figure(figsize=(4, 3))

        i = 0
        for strain in unique_strains:
            for T in unique_T:
                mask = (int_T == int(T)) & self._strain_mask(data, strain)
                n_pts = mask.sum()
                if n_pts == 0:
                    continue
                strain_str = "" if strain is None or all(s == 0 for s in strain) \
                    else f" e={strain}"
                print(f"{n_pts} data points for T = {T}{strain_str}")
                try:
                    unique_r, dF_ave, dF_var = self.compute_dF_statistics(
                        data, T, strain=strain)
                except AssertionError as e:
                    print(e)
                    continue

                dF_spline = CubicSpline(unique_r, dF_ave, bc_type="not-a-knot")
                r_dense = np.linspace(unique_r.min(), unique_r.max(), 100)
                dF_dense = dF_spline(r_dense)
                F_dense = cumulative_trapezoid(dF_dense, r_dense, initial=0)
                if symmetric:
                    F_dense += cumulative_trapezoid(dF_dense[::-1],
                                                    r_dense[::-1], initial=0)[::-1]
                    F_dense /= 2.0

                dF_var_spline = CubicSpline(unique_r, dF_var, bc_type="not-a-knot")
                F_var_dense = np.abs(
                    cumulative_trapezoid(dF_var_spline(r_dense), r_dense, initial=0)
                )
                F_var_dense += np.abs(
                    cumulative_trapezoid(dF_var_spline(r_dense)[::-1],
                                         r_dense[::-1], initial=0)
                )[::-1]
                F_var_dense /= 2.0

                label = f"F (T={T}{strain_str})"
                plt.fill_between(
                    r_dense,
                    F_dense - np.sqrt(F_var_dense),
                    F_dense + np.sqrt(F_var_dense),
                    alpha=0.2,
                    facecolor=f"C{i}",
                )
                plt.plot(r_dense, F_dense, f"C{i}-", lw=3, label=label)
                i += 1

        plt.legend()
        if filename is not None:
            plt.savefig(filename)
        else:
            plt.tight_layout()
            plt.show()

        return r_dense, F_dense

    def to_dataframe(self, data=None, symmetric=False):
        """Build a pandas DataFrame with integrated free energy per (r, T, strain).

        For each unique combination of temperature and strain, integrates
        ``dF`` over ``r`` to obtain ``F`` and propagates the block variance
        to obtain ``F_std``.

        Args:
            data: Data dictionary.  If ``None``, uses ``self.data``.
            symmetric: Average forward and backward integration.

        Returns:
            ``pandas.DataFrame`` with columns ``r``, ``T``,
            ``strain_xx``, ``strain_yy``, ``strain_zz``,
            ``strain_yz``, ``strain_xz``, ``strain_xy``,
            ``F``, ``F_std``.
        """
        import pandas as pd

        if data is None:
            data = self.data

        int_T = np.array(data["target_T"]).astype(int)
        unique_T = np.sort(np.unique(int_T))
        unique_strains = self._unique_strains(data)
        voigt_keys = ["strain_xx", "strain_yy", "strain_zz",
                      "strain_yz", "strain_xz", "strain_xy"]

        rows = []
        for strain in unique_strains:
            for T in unique_T:
                try:
                    unique_r, dF_ave, dF_var = self.compute_dF_statistics(
                        data, T, strain=strain)
                except AssertionError:
                    continue

                dF_spline = CubicSpline(unique_r, dF_ave, bc_type="not-a-knot")
                dF_var_spline = CubicSpline(unique_r, dF_var, bc_type="not-a-knot")

                F = cumulative_trapezoid(dF_spline(unique_r), unique_r, initial=0)
                if symmetric:
                    F += cumulative_trapezoid(
                        dF_spline(unique_r)[::-1], unique_r[::-1], initial=0
                    )[::-1]
                    F /= 2.0

                F_var = np.abs(
                    cumulative_trapezoid(dF_var_spline(unique_r), unique_r, initial=0)
                )
                F_var += np.abs(
                    cumulative_trapezoid(
                        dF_var_spline(unique_r)[::-1], unique_r[::-1], initial=0
                    )[::-1]
                )
                F_var /= 2.0

                s_vals = dict(zip(voigt_keys,
                                  strain if strain is not None else (0,)*6))
                for r_val, f_val, fv_val in zip(unique_r, F, F_var):
                    rows.append({"r": r_val, "T": int(T), **s_vals,
                                 "F": f_val, "F_std": np.sqrt(fv_val)})

        return pd.DataFrame(rows)


class PAFI:
    """Projected Average Force Integrator using ASE calculators.

    Builds a cubic-spline pathway from a set of NEB images and provides
    methods to run constrained Langevin MD simulations at a given reaction
    coordinate ``r`` and temperature ``T``, accumulate statistics, and
    integrate the projected force to obtain the free-energy profile.

    Attributes:
        system: Working ASE ``Atoms`` object (copy of ``images[0]``).
        calc: ASE calculator attached to ``system``.
        data: Dictionary accumulating per-step observables from :meth:`run`.
        r: Current reaction coordinate value (``0 ≤ r ≤ 1``).
        cell: Lattice matrix at the current ``r`` (shape ``(3, 3)``).
        position: Atom positions at the current ``r`` (flat array).
        tangent: Path tangent vector at the current ``r`` (flat, not normalised).
        dtangent: Second derivative of the path divided by ``|t|²`` at ``r``.
        depsilon: True strain-rate tensor ``dC/dr · C⁻¹`` at ``r``.
        E0: Spline energy at the current ``r``.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, images, calc, comm=None):
        """Initialise PAFI from a list of NEB images and an ASE calculator.

        Builds cubic splines for atomic positions, cell parameters, and
        energies from ``images``, attaches ``calc`` to the working system,
        and prepares an empty data store.

        Args:
            images: Ordered list of ASE ``Atoms`` objects spanning the
                pathway (at least two required).
            calc: ASE calculator used for forces and energies.  The same
                calculator instance is reused throughout; make sure it is
                thread-safe if you call :meth:`run` concurrently.
            comm: ASE-compatible communicator (e.g. ``MPI4PY(mpi_comm)``)
                passed to the Langevin integrator.  Controls whether random
                forces are broadcast across ranks.  ``None`` (default) uses
                ASE's default behaviour (single-process).
        """
        self.system = images[0].copy()
        self.calc = calc
        self.system.calc = self.calc
        self.comm = comm
        self.rng = np.random.default_rng()

        self.kB = units.kB
        # Inverse friction coefficient γ⁻¹; units: time / (force) = fs/(eV/Å)
        self.inv_gamma = 0.05 * units.fs / (units.eV / units.Angstrom)

        self._deformation = np.eye(3)
        self._strain_voigt = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        self._build_splines(images)
        self._data = PAFIData()

    # ------------------------------------------------------------------
    # Pathway — spline construction and navigation
    # ------------------------------------------------------------------

    def _build_splines(self, images):
        """Build cubic splines for positions, cell, energy, and strain along the path.

        Computes a dimensionless arc-length coordinate ``s ∈ [0, 1]`` from
        the Euclidean distance between PBC-wrapped displacements, then fits
        ``CubicSpline`` objects (not-a-knot boundary conditions) for:

        * atomic positions
        * cell parameters
        * NEB energies
        * strain-rate tensor ``dC/dr · C⁻¹``

        Called once during :meth:`__init__`.

        Args:
            images: Ordered list of ASE ``Atoms`` objects (≥ 2).

        Raises:
            AssertionError: If ``images`` is not a list, has fewer than two
                entries, images have different atom counts, or fewer than 5
                atoms per image.
        """
        assert isinstance(images, list), "images must be a list of ASE Atoms objects"
        assert len(images) > 1, "At least two images are required for the NEB path"
        n_atoms = [len(img) for img in images]
        assert len(set(n_atoms)) == 1, f"All images must have the same number of atoms; got {n_atoms}"
        assert n_atoms[0] >= 5, f"System must have at least 5 atoms (got {n_atoms[0]}); virial temperature requires >4 free DOF"

        self._set_pbc_function()

        spline_coord = []
        neb_energies = []
        cell_coords = []
        atom_coords = []

        for image in images:
            self.system.set_cell(image.get_cell(), scale_atoms=False)
            self.system.set_positions(image.get_positions())
            dx = self.pbc(image.get_positions() - images[0].get_positions())
            spline_coord.append(np.linalg.norm(dx))
            neb_energies.append(self.system.get_potential_energy())
            cell_coords.append(self.system.get_cell().flatten())
            atom_coords.append(images[0].get_positions().flatten() + dx)

        spline_coord = np.array(spline_coord)
        spline_coord /= spline_coord.max()

        self.spline_coordinate = spline_coord
        self.neb_energies = np.array(neb_energies)
        self.cell_coordinates = np.array(cell_coords)
        self.atom_coordinates = np.array(atom_coords)

        # Restore to initial image
        self.system.set_cell(images[0].get_cell(), scale_atoms=False)
        self.system.set_positions(images[0].get_positions())

        # Fit primary splines
        self.energy_spline = CubicSpline(
            spline_coord, self.neb_energies, bc_type="not-a-knot"
        )
        self.cell_spline = CubicSpline(
            spline_coord, self.cell_coordinates, axis=0, bc_type="not-a-knot"
        )
        self.atom_spline = CubicSpline(
            spline_coord, self.atom_coordinates, axis=0, bc_type="not-a-knot"
        )

        # True strain-rate tensor ε′ = dC/dr · C⁻¹  (shape: n_images × 3 × 3)
        dcell = self.cell_spline.derivative(nu=1)(spline_coord).reshape((-1, 3, 3))
        cells = self.cell_coordinates.reshape((-1, 3, 3))
        depsilon_coords = np.array([dc @ np.linalg.inv(c) for dc, c in zip(dcell, cells)])
        self.depsilon_spline = CubicSpline(
            spline_coord, depsilon_coords, axis=0, bc_type="not-a-knot"
        )

        self._set_pbc_function()
        self.set_r(0.0)

        self.system.set_cell(self.cell, scale_atoms=False)
        self.system.set_positions(self.position.reshape((-1, 3)))

    def _set_pbc_function(self):
        """Update ``self.pbc`` to wrap displacements into the current cell.

        ``self.pbc(v)`` maps a flat displacement vector ``v`` back into the
        primary cell of ``self.system`` using the minimum-image convention.
        Must be called whenever the cell changes.
        """
        c = self.system.get_cell()
        ic = np.linalg.inv(c)
        self.pbc = lambda v: (
            v.flatten()
            - (np.floor(np.dot(v.reshape((-1, 3)), ic) + 0.5) @ c).flatten()
        )

    def set_r(self, r):
        """Move the working system to reaction coordinate ``r``.

        Evaluates all splines at ``r`` and updates:
        ``cell``, ``dcell``, ``inv_cell``, ``V``,
        ``position``, ``N``,
        ``tangent``, ``dtangent``, ``depsilon``, ``E0``.

        Also calls ``system.set_positions`` so that the attached calculator
        sees the correct geometry.

        Args:
            r: Reaction coordinate, ``0 ≤ r ≤ 1``.

        Raises:
            AssertionError: If ``r`` is outside ``[0, 1]``.
        """
        assert 0.0 <= r <= 1.0, "r must be between 0 and 1"

        self.r = r
        F = self._deformation

        self.cell = self.cell_spline(r).reshape((3, 3)) @ F.T
        self.dcell = self.cell_spline.derivative(nu=1)(r).reshape((3, 3)) @ F.T
        self.inv_cell = np.linalg.inv(self.cell)
        self.V = np.linalg.det(self.cell)

        self.position = (self.atom_spline(r).reshape((-1, 3)) @ F.T).flatten()
        self.N = self.position.size // 3

        # Tangent: raw first derivative (not normalised; FixedMode normalises internally)
        t = (self.atom_spline.derivative(nu=1)(r).reshape((-1, 3)) @ F.T).flatten()
        norm_t = np.linalg.norm(t)
        self.tangent_sq_norm = norm_t ** 2
        self.tangent = t.copy()

        # Second derivative divided by |t|² — used in the Psi correction
        dt2 = (self.atom_spline.derivative(nu=2)(r).reshape((-1, 3)) @ F.T).flatten()
        self.dtangent = dt2 / norm_t ** 2

        # depsilon = dC/dr @ C^{-1} is invariant under uniform strain
        self.depsilon = self.depsilon_spline(r).reshape((3, 3))
        self.E0 = self.energy_spline(r)

        self.system.set_positions(self.position.reshape((-1, 3)))

    def set_strain(self, epsilon):
        """Apply a uniform strain to the pathway.

        The deformation gradient ``F = I + epsilon`` is applied to the
        cell and atomic positions when :meth:`set_r` is called.  The
        ``depsilon`` strain-rate tensor is invariant under uniform strain.

        Args:
            epsilon: 3x3 strain tensor (symmetric, small-strain convention).
                Use ``np.zeros((3, 3))`` to reset.
        """
        epsilon = np.asarray(epsilon, dtype=float)
        assert epsilon.shape == (3, 3), "epsilon must be a 3x3 array"
        self._deformation = np.eye(3) + epsilon
        # Store Voigt notation (xx, yy, zz, yz, xz, xy) as hashable tuple
        self._strain_voigt = (
            epsilon[0, 0], epsilon[1, 1], epsilon[2, 2],
            epsilon[1, 2], epsilon[0, 2], epsilon[0, 1],
        )

    # ------------------------------------------------------------------
    # PAFI force projection
    # ------------------------------------------------------------------

    def pafi_projection(self, atoms, T):
        """Compute PAFI observables for ``atoms`` at temperature ``T``.

        Projects the atomic forces and virial stress onto the path tangent and
        applies the Psi curvature correction to obtain the instantaneous
        estimate of ``dF/dr``.

        Args:
            atoms: ASE ``Atoms`` object with a live calculator attached.
            T: Thermostat temperature in Kelvin (used for the stress
                entropic correction).

        Returns:
            Tuple ``(dx, dF, Psi, virialT, dxt, dt)`` where:

            * ``dx`` — PBC-wrapped displacement from path position (flat array, Å).
            * ``dF`` — Instantaneous projected free-energy gradient (eV).
            * ``Psi`` — Curvature correction factor (dimensionless).
            * ``virialT`` — Virial-theorem kinetic temperature estimate (K).
            * ``dxt`` — Component of ``dx`` along the tangent (Å).
            * ``dt`` — Wall-clock time for this call (s).
        """
        t0 = time.time()

        x = atoms.get_positions().flatten()
        dx = self.pbc(x - self.position.flatten())
        f = self.calc.get_forces(atoms).flatten()

        # Reconstruct 3×3 stress tensor from Voigt notation [xx,yy,zz,yz,xz,xy]
        _s = self.calc.get_stress(atoms)
        s = np.zeros((3, 3))
        s[0, 0], s[1, 1], s[2, 2] = _s[0], _s[1], _s[2]
        s[0, 1] = s[1, 0] = _s[3]
        s[0, 2] = s[2, 0] = _s[4]
        s[1, 2] = s[2, 1] = _s[5]

        # Projections along the path tangent
        dx_path = np.dot(dx, self.tangent)
        f_path = np.dot(f, self.tangent)

        # Stress contribution (mechanical + entropic volume term)
        s_path = -np.dot(s.flatten(), self.depsilon.flatten()) * atoms.get_volume()
        s_path -= self.kB * T * (self.dcell.flatten() @ self.inv_cell.flatten()) * self.N

        # Psi: curvature correction (≈1 near the path; < 1 in high-curvature regions)
        Psi = 1.0 - np.dot(dx, self.dtangent)
        dF = Psi * (f_path + s_path)

        # Virial temperature (excluding ~4 constrained DOF); undefined for ≤4 atoms
        n_free = 3.0*len(self.system) - 4.0
        virialT = -np.dot(dx, f) / self.kB / n_free if n_free > 0 else np.nan

        dxt = np.dot(dx, self.tangent)
        dt = time.time() - t0

        return dx, dF, Psi, virialT, dxt, dt

    # ------------------------------------------------------------------
    # Simulation
    # ------------------------------------------------------------------

    def run(
        self,
        r=0.0,
        T=300.0,
        nsteps=1000,
        thermsteps=500,
        minsteps=5,
        verbose=False,
        print_nevery=50,
        friction=0.1,
        postmin=False,
        return_dx=False,
    ):
        """Run a constrained Langevin simulation at reaction coordinate ``r``.

        Thermalises the system for ``thermsteps`` steps, then collects PAFI
        observables every MD step for ``nsteps`` steps.  Results are appended
        to ``self.data``.  Each call receives a UUID barcode so that
        independent runs at the same ``(r, T)`` can be identified.

        Args:
            r: Reaction coordinate value, ``0 ≤ r ≤ 1``.  Default ``0.0``.
            T: Target temperature in Kelvin.  ``T=0`` reduces to a single
                force evaluation with no thermalisation.  Default ``300.0``.
            nsteps: Number of production MD steps.  Default ``1000``.
            thermsteps: Number of Langevin thermalisation steps before data
                collection.  Default ``500``.
            minsteps: Maximum number of FIRE2 minimisation steps applied
                after the run (only when ``postmin=True``).  Default ``5``.
            verbose: Print a header and per-step summary if ``True``.
                Default ``False``.
            print_nevery: Print frequency when ``verbose=True``.
                Default ``5``.
            friction: Langevin friction coefficient (ASE units).
                Default ``0.1``.
            postmin: Run a short FIRE2 relaxation of ``self.system`` after
                data collection.  Useful for seeding the next :meth:`set_r`
                call.  Default ``True``.
            return_dx: If ``True``, return the mean-subtracted displacement
                field ``dx`` (shape ``(N, 3)``) of the post-minimisation
                system.  Default ``False``.

        Returns:
            If ``return_dx`` is ``True``, returns a ``(N, 3)`` NumPy array
            of the displacement from the path reference position after
            minimisation.  Otherwise returns ``None``.

        Raises:
            AssertionError: If ``nsteps ≤ 0``, ``thermsteps ≤ 0``, or
                ``T < 0``.
        """
        assert nsteps > 0, "nsteps must be greater than 0"
        assert thermsteps > 0, "thermsteps must be greater than 0"
        assert T >= 0.0, "Temperature must be non-negative"

        pseudo_unique_hex = uuid.uuid4().int

        if T == 0.0:
            nsteps = 1
            thermsteps = 0

        self.set_r(r)
        self.system.set_cell(self.cell, scale_atoms=False)
        self.system.set_positions(self.position.reshape((-1, 3)))

        atoms = self.system.copy()
        atoms.calc = self.calc
        atoms.set_constraint(FixedMode(self.tangent))

        langevin_kw = dict(atoms=atoms, timestep=units.fs, temperature_K=T,
                           friction=friction, rng=self.rng)
        if self.comm is not None:
            langevin_kw['comm'] = self.comm
        dyn = Langevin(**langevin_kw)
        t0 = time.time()
        dyn.run(thermsteps)
        if verbose:
            print(f"Thermalisation over {thermsteps} steps took {time.time()-t0:.2f}s")
            print(
                f"{'Step':<8}{'r':<8}{'T':<8}{'dF':<8}"
                f"{'Psi':<8}{'dx.t':<8}{'max(dx)':<8}{'dt':<6}"
            )

        # Production loop
        for step in range(nsteps):
            dyn.run(1)
            dx, dF, Psi, virialT, dxt, dt = self.pafi_projection(atoms, T)
            maxdx = float(np.max(np.abs(dx)))

            self.data["barcode"].append(pseudo_unique_hex)
            self.data["step"].append(step)
            self.data["r"].append(self.r)
            self.data["T"].append(virialT)
            self.data["target_T"].append(T)
            self.data["strain"].append(list(self._strain_voigt))
            self.data["dF"].append(float(dF))
            self.data["Psi"].append(float(Psi))
            self.data["dx.t"].append(float(dxt))
            self.data["max(dx)"].append(maxdx)
            self.data["dt"].append(dt)

            if verbose and step % print_nevery == 0:
                print(
                    f"{step:<8}{self.r:<8.2f}{virialT:<8.2f}{dF:<8.2f}"
                    f"{Psi:<8.2f}{dxt:<8.2f}{maxdx:<8.2f}{dt:<6.3f}s"
                )

        if postmin:
            opt = FIRE2(atoms=self.system)
            opt.run(fmax=0.015, steps=minsteps)

        if return_dx:
            dx = self.pbc(
                self.system.get_positions().flatten() - self.position
            ).reshape((-1, 3))
            dx -= dx.mean(0)
            return dx

    # ------------------------------------------------------------------
    # Data management (delegated to PAFIData)
    # ------------------------------------------------------------------

    @property
    def data(self):
        """Data dictionary accumulating per-step observables."""
        return self._data.data

    @data.setter
    def data(self, value):
        self._data.data = value

    def reset_data(self):
        """Clear all accumulated simulation data."""
        self._data.reset_data()

    def write_json(self, filename):
        """Serialise ``self.data`` to a JSON file."""
        self._data.write_json(filename)

    def read_json(self, filename):
        """Load simulation data from a JSON file."""
        self._data.read_json(filename)

    def compute_dF_statistics(self, data, T, strain=None):
        """Compute mean and variance of the projected force at each ``r``."""
        return self._data.compute_dF_statistics(data, T, strain=strain)

    def plot_data(self, data=None, filename=None, symmetric=False):
        """Plot the free-energy profile ``F(r)``."""
        return self._data.plot_data(data=data, filename=filename, symmetric=symmetric)

    # ------------------------------------------------------------------
    # Plotting (requires calculator and splines)
    # ------------------------------------------------------------------

    def plot_neb(self, filename=None):
        """Plot the NEB reference path and the PAFI force-work integral.

        Evaluates the spline energies and projected forces at each knot, fits
        a work spline ``W(r) = -∫ F·dX₀``, and overlays:

        * knot energies (scatter),
        * energy spline (dashed),
        * work integral (line + markers).

        Args:
            filename: If given, save the figure to this path instead of
                displaying it interactively.
        """
        E_spl, dE_spl = [], []
        for s in self.spline_coordinate:
            self.set_r(s)
            self.system.set_cell(self.cell, scale_atoms=False)
            self.system.set_positions(self.position.reshape((-1, 3)))
            dF = -self.pafi_projection(self.system, 0.0)[1]
            E_spl.append(self.system.get_potential_energy())
            dE_spl.append(dF)

        self.set_r(0.0)

        E_spl = np.array(E_spl)
        dE_spl = np.array(dE_spl)
        work_spline = CubicSpline(self.spline_coordinate, dE_spl, bc_type="not-a-knot")
        r = np.linspace(self.spline_coordinate.min(), self.spline_coordinate.max(),
                        2 * len(self.spline_coordinate))
        work = cumulative_trapezoid(work_spline(r), r, initial=0)

        plt.figure(figsize=(4, 3))
        plt.plot(self.spline_coordinate, E_spl - E_spl[0], "o", label="Knot Energies")
        plt.plot(r, self.energy_spline(r) - E_spl[0], "--", label="Energy Spline")
        plt.plot(r, work, ".-",
                 label=r"$-\int_0^r{\bf F}\cdot{\rm d}{\bf X}_0(r)$ (Work)")
        plt.legend()

        if filename is not None:
            plt.savefig(filename)


class PAFIMPI:
    """PAFI sampling distributed over MPI workers.

    Each rank owns an independent :class:`PAFI` instance with its own
    calculator.  All ranks call :meth:`run` (or :meth:`run_distributed`)
    collectively; the resulting data are gathered to rank 0 of
    ``inter_comm``, which then holds the full combined sample for writing
    and analysis.

    Args:
        images: Ordered list of ASE ``Atoms`` objects spanning the pathway
            (passed unchanged to :class:`PAFI` on every rank).
        calc: ASE calculator **or** list of ASE calculators.
            * Single calculator — every rank uses the same (independently
              initialised) calculator instance, which is the normal case when
              each MPI process already owns its own object (e.g. one
              ``LAMMPSlib`` or ``MACECalculator`` per rank).
            * List of calculators — rank *i* uses ``calc[i % len(calc)]``,
              which is convenient when you have fewer devices than ranks
              (e.g. two GPUs shared across eight processes).

            **Important:** the calculator must be constructed with the same
            communicator as ``worker_comm``.  There is no generic ASE
            interface for extracting or overriding a calculator's
            communicator, so this is the caller's responsibility.  For
            example, ``LAMMPSlib(comm=worker_comm, ...)`` or
            ``MACECalculator(comm=worker_comm, ...)``.
        worker_comm: MPI communicator for each worker's internal parallelism.
            Passed to the Langevin integrator so that each worker group
            generates independent random forces.  Defaults to
            ``MPI.COMM_SELF`` (one process per worker).  Set to a
            sub-communicator for multi-process workers (e.g. 4 CPU cores per
            LAMMPS instance, or 2 GPUs per worker).
        inter_comm: MPI communicator across workers, used for gather/scatter
            of PAFI data.  Only the root rank of each ``worker_comm`` should
            participate. Defaults to ``MPI.COMM_WORLD``.

    Example — 1 process per worker (default)::

        pafi = PAFIMPI(images, calc)

    Example — 4 CPU cores per worker::

        worker_comm = MPI.COMM_WORLD.Split(rank // 4, rank)
        is_root = (rank % 4 == 0)
        inter_comm = MPI.COMM_WORLD.Split(0 if is_root else MPI.UNDEFINED, rank // 4)
        pafi = PAFIMPI(images, calc, worker_comm=worker_comm, inter_comm=inter_comm)

    Example — all ranks sample the same (r, T) point::

        pafi = PAFIMPI(images, calc)
        for T in temperatures:
            for r in np.linspace(0, 1, 11):
                pafi.run(r, T, nsteps=1000, thermsteps=500)
                pafi.write_json(f"T{int(T)}_r{int(100*r)}.json")
                pafi.reset_data()

    Example — distribute r values across ranks::

        pafi = PAFIMPI(images, calc)
        for T in temperatures:
            pafi.run_distributed(np.linspace(0, 1, 11), T, nsteps=1000)
            pafi.write_json(f"T{int(T)}_all_r.json")
            pafi.reset_data()
    """

    def __init__(self, images, calc, worker_comm=None, inter_comm=None):
        from mpi4py import MPI
        from ase.parallel import MPI4PY

        if worker_comm is None:
            worker_comm = MPI.COMM_SELF
        if inter_comm is None:
            inter_comm = MPI.COMM_WORLD

        self.worker_comm = worker_comm
        self.is_worker_root = (worker_comm.Get_rank() == 0)

        # inter_comm may be MPI.COMM_NULL for non-root worker ranks
        if self.is_worker_root:
            self.inter_comm = inter_comm
            self.rank = inter_comm.Get_rank()
            self.size = inter_comm.Get_size()
        else:
            self.inter_comm = None
            self.rank = None
            self.size = None

        if isinstance(calc, list):
            _calc = calc[self.rank % len(calc)] if self.is_worker_root else calc[0]
        else:
            _calc = calc

        self._pafi = PAFI(images, _calc, comm=MPI4PY(worker_comm))
        self._history = PAFIData()

    # ------------------------------------------------------------------
    # Collective sampling
    # ------------------------------------------------------------------

    @property
    def data(self):
        """Full data on rank 0 (history + last gather); local data on others."""
        return self._history.data

    def run(
        self,
        r,
        T,
        nsteps=1000,
        thermsteps=500,
        minsteps=5,
        verbose=False,
        print_nevery=50,
        friction=0.1,
        postmin=False,
        return_dx=False,
        gather=True,
    ):
        """Run PAFI on every rank at the same ``(r, T)`` point.

        Each rank performs an independent Langevin trajectory.  When
        ``gather=True`` (default) the per-rank data are merged on rank 0 so
        that :meth:`write_json` immediately writes the combined sample.

        Args:
            r, T, nsteps, thermsteps, minsteps, verbose, print_nevery,
            friction, postmin, return_dx: forwarded to :meth:`PAFI.run`.
            gather: If ``True``, call :meth:`_gather_data` after the run.

        Returns:
            The return value of the underlying :meth:`PAFI.run` call on this
            rank (only meaningful when ``return_dx=True``).
        """
        result = self._pafi.run(
            r,
            T,
            nsteps=nsteps,
            thermsteps=thermsteps,
            minsteps=minsteps,
            verbose=verbose,
            print_nevery=print_nevery,
            friction=friction,
            postmin=postmin,
            return_dx=return_dx,
        )
        if gather:
            self._gather_data()
        return result

    def run_distributed(
        self,
        r_values,
        T,
        nsteps=1000,
        thermsteps=500,
        minsteps=5,
        verbose=False,
        print_nevery=50,
        friction=0.1,
        postmin=False,
        gather=True,
    ):
        """Distribute ``r_values`` across ranks and run PAFI at each assigned point.

        Rank *i* processes ``r_values[i::size]`` (stride distribution, so the
        load is spread evenly and the r values remain well-interleaved across
        ranks).  When ``gather=True`` (default) all per-rank data are gathered
        to rank 0 after the local runs complete.

        Args:
            r_values: Sequence of reaction-coordinate values to sample.
            T: Temperature in K (scalar, applied to every r value).
            nsteps, thermsteps, minsteps, verbose, print_nevery, friction,
            postmin: forwarded to :meth:`PAFI.run`.
            gather: If ``True``, call :meth:`_gather_data` after all local runs.
        """
        # Worker root determines its assigned r values; broadcast to all
        # ranks in worker_comm so non-roots follow in lockstep (required for
        # MPI-parallel calculators where get_forces is collective).
        if self.is_worker_root:
            my_r = list(r_values[self.rank :: self.size])
        else:
            my_r = None
        my_r = self.worker_comm.bcast(my_r, root=0)
        for r in my_r:
            self._pafi.run(
                r,
                T,
                nsteps=nsteps,
                thermsteps=thermsteps,
                minsteps=minsteps,
                verbose=verbose,
                print_nevery=print_nevery,
                friction=friction,
                postmin=postmin,
            )
        if gather:
            self._gather_data()

    # ------------------------------------------------------------------
    # Data management
    # ------------------------------------------------------------------

    def _gather_data(self):
        """Gather fresh per-rank data and append to rank 0's history.

        Each rank's ``_pafi.data`` (samples produced since the last gather)
        is collected via ``inter_comm``.  On rank 0 the gathered samples are
        appended to ``_history``.  All ranks then reset ``_pafi.data`` so
        the next run starts with a clean buffer.
        """
        if not self.is_worker_root:
            self._pafi.reset_data()
            return
        all_data = self.inter_comm.gather(self._pafi.data, root=0)
        if self.rank == 0:
            for rank_data in all_data:
                for key in self._history.data:
                    self._history.data[key].extend(rank_data[key])
        self._pafi.reset_data()

    def reset_data(self):
        """Reset accumulated data (history and pending) on all ranks."""
        self._pafi.reset_data()
        self._history.reset_data()

    def write_json(self, filename):
        """Write gathered data to JSON (inter_comm rank 0 only)."""
        if self.is_worker_root and self.rank == 0:
            self._history.write_json(filename)

    def read_json(self, filename):
        """Read JSON into rank 0's history."""
        if self.is_worker_root and self.rank == 0:
            self._history.read_json(filename)

    # ------------------------------------------------------------------
    # Analysis and plotting (inter_comm rank 0 delegates to history)
    # ------------------------------------------------------------------

    def compute_dF_statistics(self, data=None, T=None, strain=None):
        """Compute dF statistics from gathered history."""
        return self._history.compute_dF_statistics(
            data=data if data is not None else self._history.data,
            T=T, strain=strain)

    def plot_data(self, data=None, filename=None):
        """Plot free-energy profile (inter_comm rank 0 only)."""
        if self.is_worker_root and self.rank == 0:
            self._history.plot_data(data=data, filename=filename)

    def plot_neb(self, filename=None):
        """Plot NEB reference path (inter_comm rank 0 only)."""
        if self.is_worker_root and self.rank == 0:
            self._pafi.plot_neb(filename=filename)
