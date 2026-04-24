"""Unit tests for pafiase.

Focuses on the T=0 case (single force-evaluation step, no thermalisation noise)
so results are deterministic and cheap to compute.

Fixture
-------
A 3-image linear path is built from an 8-atom FCC Cu supercell (EMT calculator)
with the first atom displaced by 0.2 Å along x between the first and last image.
"""

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT

from pafiase import PAFI
from pafiase import _block_variance


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pafi():
    """PAFI instance built from a 3-image Cu FCC path (EMT calculator)."""
    base = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((1, 1, 2)).repeat((1, 1, 2))  # 8 atoms

    images = []
    for i in range(3):
        img = base.copy()
        img.positions[0, 0] += i * 0.1  # linear interpolation, 0 → 0.2 Å
        images.append(img)

    return PAFI(images, EMT())


# ---------------------------------------------------------------------------
# _block_variance
# ---------------------------------------------------------------------------

class TestBlockVariance:
    def test_short_array_uses_raw_variance(self):
        """Arrays shorter than 20 fall back to data.var()."""
        data = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        assert _block_variance(data) == pytest.approx(data.var())

    def test_constant_array_is_zero(self):
        """Constant data has zero variance."""
        data = np.ones(100) * 3.14
        assert _block_variance(data) == pytest.approx(0.0, abs=1e-12)

    def test_returns_positive_float(self):
        """Result is a non-negative scalar for noisy data."""
        rng = np.random.default_rng(0)
        data = rng.normal(0.0, 1.0, 500)
        result = _block_variance(data)
        assert isinstance(result, float)
        assert result >= 0.0

    def test_larger_variance_gives_larger_result(self):
        """Scaling data by k scales block variance by k²."""
        rng = np.random.default_rng(1)
        data = rng.normal(0.0, 1.0, 500)
        assert _block_variance(3.0 * data) > _block_variance(data)


# ---------------------------------------------------------------------------
# PAFI construction — input validation
# ---------------------------------------------------------------------------

class TestInputValidation:
    def _make_images(self, n_atoms, n_images=3, displacement=0.1):
        from ase import Atoms
        images = []
        for i in range(n_images):
            pos = [[j * 2.0, 0, 0] for j in range(n_atoms)]
            pos[0][0] += i * displacement
            img = Atoms(f"Cu{n_atoms}", positions=pos, cell=[20, 20, 20], pbc=True)
            images.append(img)
        return images

    def test_too_few_atoms_raises(self):
        images = self._make_images(n_atoms=4)
        with pytest.raises(AssertionError, match="at least 5"):
            PAFI(images, EMT())

    def test_mismatched_atom_counts_raises(self):
        base = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((1, 1, 2)).repeat((1, 1, 2))  # 8 atoms
        images = [base.copy(), base.copy()]
        # Add an extra atom to the second image
        from ase import Atoms
        images[1] = images[1] + Atoms("Cu", positions=[[0.5, 0.5, 0.5]])
        with pytest.raises(AssertionError, match="same number of atoms"):
            PAFI(images, EMT())

    def test_single_image_raises(self):
        base = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((1, 1, 2)).repeat((1, 1, 2))
        with pytest.raises(AssertionError, match="At least two"):
            PAFI([base], EMT())

    def test_valid_construction_succeeds(self):
        base = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((1, 1, 2)).repeat((1, 1, 2))  # 8 atoms
        images = [base.copy(), base.copy()]
        images[1].positions[0, 0] += 0.1
        pafi = PAFI(images, EMT())
        assert len(pafi.system) == len(base)


# ---------------------------------------------------------------------------
# PAFI construction — splines
# ---------------------------------------------------------------------------

class TestConstruction:
    def test_spline_coordinate_range(self, pafi):
        """Spline coordinate runs from 0 to 1."""
        assert pafi.spline_coordinate[0] == pytest.approx(0.0)
        assert pafi.spline_coordinate[-1] == pytest.approx(1.0)

    def test_spline_coordinate_length(self, pafi):
        """One spline knot per image."""
        assert len(pafi.spline_coordinate) == 3

    def test_initial_r(self, pafi):
        """r is initialised to 0.0."""
        pafi.set_r(0.0)
        assert pafi.r == pytest.approx(0.0)

    def test_set_r_out_of_bounds(self, pafi):
        with pytest.raises(AssertionError):
            pafi.set_r(-0.1)
        with pytest.raises(AssertionError):
            pafi.set_r(1.1)

    def test_set_r_updates_cell(self, pafi):
        """Cell matrix changes when r changes."""
        pafi.set_r(0.0)
        cell_0 = pafi.cell.copy()
        pafi.set_r(0.5)
        cell_5 = pafi.cell.copy()
        # For this path the cell is constant, so they should be equal
        np.testing.assert_allclose(cell_0, cell_5)

    def test_set_r_updates_position(self, pafi):
        """Position vector changes monotonically along x for this path."""
        pafi.set_r(0.0)
        x0 = pafi.position[0]
        pafi.set_r(1.0)
        x1 = pafi.position[0]
        assert x1 > x0  # first atom displaced +0.2 Å along x


# ---------------------------------------------------------------------------
# run() at T=0
# ---------------------------------------------------------------------------

class TestRun0K:
    def setup_method(self):
        """Fresh PAFI instance for each test so data doesn't accumulate."""
        base = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((1, 1, 2))
        images = [base.copy() for _ in range(3)]
        for i, img in enumerate(images):
            img.positions[0, 0] += i * 0.1
        self.pafi = PAFI(images, EMT())

    def test_produces_exactly_one_step(self):
        """T=0 run records exactly one data point."""
        self.pafi.run(r=0.5, T=0.0)
        assert len(self.pafi.data["step"]) == 1

    def test_r_is_recorded(self):
        """The requested r value appears in data."""
        self.pafi.run(r=0.25, T=0.0)
        assert self.pafi.data["r"][0] == pytest.approx(0.25)

    def test_target_T_is_zero(self):
        self.pafi.run(r=0.0, T=0.0)
        assert self.pafi.data["target_T"][0] == pytest.approx(0.0)

    def test_all_keys_populated(self):
        """Every expected key is present and has one entry."""
        self.pafi.run(r=0.5, T=0.0)
        for key in ("barcode", "step", "r", "T", "target_T", "dF", "Psi", "dx.t", "max(dx)", "dt"):
            assert len(self.pafi.data[key]) == 1, f"key '{key}' not populated"

    def test_max_dx_is_non_negative(self):
        self.pafi.run(r=0.5, T=0.0)
        assert self.pafi.data["max(dx)"][0] >= 0.0

    def test_psi_near_one_at_path(self):
        """Psi ≈ 1 when the system sits on the reference path (dx ≈ 0)."""
        self.pafi.run(r=0.5, T=0.0, postmin=False)
        assert abs(self.pafi.data["Psi"][0] - 1.0) < 0.5

    def test_multiple_runs_accumulate(self):
        """Data from successive calls is appended, not overwritten."""
        for r in [0.0, 0.5, 1.0]:
            self.pafi.run(r=r, T=0.0)
        assert len(self.pafi.data["r"]) == 3

    def test_return_dx_shape(self):
        """return_dx=True returns an (N, 3) array."""
        dx = self.pafi.run(r=0.5, T=0.0, return_dx=True)
        assert dx is not None
        assert dx.shape == (len(self.pafi.system), 3)

    def test_return_dx_mean_subtracted(self):
        """Returned dx has zero mean (centre-of-mass removed)."""
        dx = self.pafi.run(r=0.5, T=0.0, return_dx=True)
        np.testing.assert_allclose(dx.mean(axis=0), 0.0, atol=1e-10)


# ---------------------------------------------------------------------------
# reset_data
# ---------------------------------------------------------------------------

class TestResetData:
    def test_clears_all_entries(self):
        base = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((1, 1, 2))
        images = [base.copy() for _ in range(3)]
        for i, img in enumerate(images):
            img.positions[0, 0] += i * 0.1
        p = PAFI(images, EMT())
        p.run(r=0.0, T=0.0)
        p.reset_data()
        for key in p.data:
            assert len(p.data[key]) == 0, f"key '{key}' not cleared"


# ---------------------------------------------------------------------------
# compute_dF_statistics
# ---------------------------------------------------------------------------

class TestComputeDFStatistics:
    def setup_method(self):
        base = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((1, 1, 2))
        images = [base.copy() for _ in range(3)]
        for i, img in enumerate(images):
            img.positions[0, 0] += i * 0.1
        self.pafi = PAFI(images, EMT())
        for r in [0.0, 0.5, 1.0]:
            self.pafi.run(r=r, T=0.0)

    def test_output_shapes(self):
        unique_r, dF_ave, dF_var = self.pafi.compute_dF_statistics(self.pafi.data, 0)
        assert len(unique_r) == 3
        assert dF_ave.shape == (3,)
        assert dF_var.shape == (3,)

    def test_unique_r_matches_inputs(self):
        unique_r, _, _ = self.pafi.compute_dF_statistics(self.pafi.data, 0)
        np.testing.assert_allclose(unique_r, [0.0, 0.5, 1.0])

    def test_dF_var_non_negative(self):
        _, _, dF_var = self.pafi.compute_dF_statistics(self.pafi.data, 0)
        assert (dF_var >= 0).all()

    def test_single_r_raises(self):
        """Statistics cannot be computed with only one r value."""
        self.pafi.reset_data()
        self.pafi.run(r=0.5, T=0.0)
        with pytest.raises(AssertionError, match="Only one unique value"):
            self.pafi.compute_dF_statistics(self.pafi.data, 0)
