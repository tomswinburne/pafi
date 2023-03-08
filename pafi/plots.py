import numpy as np
from glob import glob
import pafi
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

def plot_neb_results(args):
    """
    Utility to plot the results of a NEB calculation from a log.lammps file.
    """
    if args.terminal:
        import plotext as plx

    d = np.genfromtxt(args.log_lammps, skip_header=3, invalid_raise=False)

    if len(np.where(np.isnan(d[:, 0]))[0])>0:
        # seems needed for some LAMMPS versions
        i_before_climb = np.where(np.isnan(d[:, 0]))[0][0]-1
    else: 
        i_before_climb = None

    d = d[np.isfinite(d).any(axis=1)]
    a = d[:, 9:].reshape((d.shape[0], int((d.shape[1]-9)/2), 2))
    if int(args.cut_at_image) > -1:
        a = a[:, :int(args.cut_at_image), :]

    cols = 1
    if args.energy_diff: 
        cols += 1
    if args.rcdiff:
        cols += 1
    
    fig, ax = plt.subplots(1, cols, figsize=(cols*int(args.size), int(args.size)))
    if cols == 1: ax = [ax]
    E_ref = a[-1, 0, 1]

    for i, (p, l) in enumerate(zip([0, i_before_climb, -1], ("Initial", "Before Climb", "Final "))):
        if p is None: continue
        if (args.final) and ("Final" not in l): continue

        norm = a[p, :, 0].max() if args.normalize_r else 1 
        n_img = a[p, :, 0].shape[0]

        if args.image_number:
            r = np.arange(len(a[p, :,1]))+1
            for axi in ax:
                axi.set_xlim(xmin=r[0]-0.2, xmax=r[-1]+0.2)
        else:
            r = a[p, :, 0]/norm

        # 1st panel: energy profile
        l += f" max. {(a[p, :, 1]-E_ref).max():.3f} eV"
        ax[0].plot(r, a[p, :, 1]-E_ref, ".-", label=l, color=f"C{i}")
        if "Final" in l:
            ax[0].fill_between(r, (a[p, :, 1]-E_ref).min(), a[p, :, 1]-E_ref, color=f"C{i}", alpha=0.1)
        # ax[0].annotate(f"{l} max: {(a[p, :, 1]-E_ref).max():.3f} eV", (5, 200-12*i), xycoords="axes points")
        if args.energy_diff:
            # 2nd panel: energy spacing
            ax[1].step(np.arange(n_img)+1, 
                    np.abs(np.diff(a[p, :, 1], prepend=a[p, 0, 1])),  color=f"C{i}")
            ax[1].fill_between(np.arange(n_img)+1, 
                    0, np.abs(np.diff(a[p, :, 1], prepend=a[p, 0, 1])), 
                    label=l, color=f"C{i}", step="pre", alpha=0.1)

            ax[1].set_title("Energy spacing")    
            ax[1].set_xlabel("Image Number")
            ax[1].set_ylabel("Energy difference (eV)")

        if args.print_profile:
            print(l, "\n", f"{a[p, :, 0]/norm}", a[p, :, 1]-a[p, 0, 1])
        if args.rcdiff:
            # 3rd panel with rc spacing
            ax[-1].step(np.arange(n_img)+1, np.abs(np.diff(a[p, :, 0], prepend=np.diff(a[p, :, 0])[0])), color=f"C{i}")
            ax[-1].fill_between(np.arange(n_img)+1, 0, np.abs(np.diff(a[p, :, 0], prepend=np.diff(a[p, :, 0])[0])), 
                                step="pre", color=f"C{i}", alpha=0.1)
            ax[-1].set_title("Reaction Coordinate spacing")
            ax[-1].set_xlabel("Image Number")
            ax[-1].set_ylabel("Reaction coordinate")

    ax[0].set_title("Energy profile")
    if args.image_number:
        ax[0].set_xlabel("Image Number")
    else:
        ax[0].set_xlabel("Reaction Coordinate")
    ax[0].set_ylabel("Energy (eV)")
    ax[0].legend()

    fig.tight_layout()
    if args.save is not None:
        try:
            plt.savefig(args.save)
        except:
            print(f"failed saving plot to {args.f}")
    if not args.quiet:
        if args.terminal:
            plx.from_matplotlib(fig)
            plx.show()
        else:
            plt.show()


def guess_mode(x):
    if (not isinstance(x, str)) or ("*" in x):
        return "m"
    else :
        return "s"


def plot_profile(args):

    cols = 1
    if args.hist:
        cols += 1
    if args.std:
        cols += 1

    fig,axs = plt.subplots(1, cols, figsize=(cols*args.size, args.size))
    if cols == 1: axs = [axs]

    plt.sca(axs[0])
    for file in glob(args.file):
        r = pafi.PafiResult(file)
        r.plot(ax=axs[0])
    if cols>1:
        try:
            # will fail at 0K
            if args.hist:
                plt.sca(axs[1])
                plt.hist(r.hist_data, bins=30)
                plt.ylabel("(dF-dFmean)/std")

            if args.std:
                plt.sca(axs[-1])
                plt.plot(r.std, "o-")
                plt.ylabel("std")
                plt.xlabel("hyperplane")
        except: pass
    plt.tight_layout()
    plt.show()



def plot_pafi_results(args):

    if args.terminal:
        import plotext as plx
    
    
    if args.mode is None:
        # If not specified, try to guess from pattern
        args.mode = guess_mode(args.file)

    # Plot a single profile 
    if args.mode == "s":
        plot_profile(args)       
    elif args.mode == "m":
        fig, ax = plt.subplots(1, 2, figsize=(args.size*2, args.size*1),dpi=144,sharey=True)

        if args.H0 is not None:
            additional_pts = [(0, args.H0)]
        else: additional_pts = None

        pafi.free_energy_vs_temperature(glob(args.file), 
                                        ax=ax,
                                        fit_harmonic=args.no_fit,
                                        add_pts=additional_pts,
                                        harmonic_until_T=args.harmonic_limit)

    fig.tight_layout()
    if args.save is not None:
        try:
            plt.savefig(args.save)
        except:
            print(f"failed saving plot to {args.f}")
    if not args.quiet:
        if args.terminal:
            plx.from_matplotlib(fig)
            plx.show()
        else:
            plt.show()


        