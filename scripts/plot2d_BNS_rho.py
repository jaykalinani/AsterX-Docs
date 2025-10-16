#!/usr/bin/env python3
import os, re, glob, shutil
import numpy as np
import imageio.v2 as imageio
import openpmd_api as io

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm
from matplotlib import cm

# ------------------ styling: LaTeX only if available ------------------
USE_TEX = shutil.which("latex") is not None
plt.rcParams.update({
    "text.usetex": bool(USE_TEX),
    "font.family": "serif" if USE_TEX else "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 10,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.labelpad": 0.5,       # label ↔ tick-label gap
    "ytick.major.pad": 0.5,     # tick-label ↔ axis spine gap
    "ytick.minor.pad": 0.5,
})

# ------------------ config ------------------
DATA_DIR   = "/scratch/09228/jkalinan/simulations/BNS_IG_fixedGrid/BNS_IG_fixedGrid"
OUT_DIR    = "/work/09228/jkalinan/vista/scripts/plots/movies/BNS_rho_frames"
MOVIE_FILE = "/work/09228/jkalinan/vista/scripts/plots/movies/BNS_rho.mp4"
os.makedirs(OUT_DIR, exist_ok=True)

# ---- NEW: merger timing controls ----
MERGER_TIME_MS = 15.0   # shift printed time: t_display = t_ms - MERGER_TIME_MS
FINAL_AFTER_MS = 30.0   # stop when t_display exceeds this (set None to disable)

# units & scales
CU_TO_KM = 1.4767161818921162  # km per CU
RHO_CU_TO_CGS = 6.1762691458861658e17     # g/cm^3 per CU
TIME_CU_TO_MS = 4.9257949707731345e-03    # ms per CU
RHO_VMIN, RHO_VMAX = (1e-11*RHO_CU_TO_CGS, 1e-3*RHO_CU_TO_CGS)

# original CU domain
XMIN_CU, XMAX_CU = -50.0, 50.0
YMIN_CU, YMAX_CU = -50.0, 50.0
ZMIN_CU, ZMAX_CU = -50.0, 50.0

# plot in kilometers
XMIN, XMAX = XMIN_CU*CU_TO_KM, XMAX_CU*CU_TO_KM
YMIN, YMAX = YMIN_CU*CU_TO_KM, YMAX_CU*CU_TO_KM
ZMIN, ZMAX = ZMIN_CU*CU_TO_KM, ZMAX_CU*CU_TO_KM

EXTENT_XY = (XMIN, XMAX, YMIN, YMAX)
EXTENT_XZ = (XMIN, XMAX, ZMIN, ZMAX)

FPS   = 12
NXNY  = 1024
CMAP = cm.get_cmap("plasma").copy()
CMAP.set_bad(color="#0b0e2c", alpha=1.0)   # NaNs -> dark ink
CMAP.set_under(color="#0b0e2c", alpha=1.0) # <vmin -> dark too

REC_COMP    = "hydrobasex_rho"
REC_NAME_RE = re.compile(r"^hydrobasex_rho_patch(\d+)_lev(\d+)$")
EDGE_FILL_PIX = 3  # try 2–5; increase if you still see a rim

try:
    from scipy.interpolate import RegularGridInterpolator
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

# ---- NEW: how aggressively to hide the refinement edge (pixels) ----
EDGE_ERODE = 1  # 1 is usually enough; increase to 2 if you still see an outline

# ------------------ helpers ------------------
def _clean_var(s):
    # strip any leading/trailing $ so we don't end up with $$y$$
    import re as _re
    return _re.sub(r'^\$+|\$+$', '', str(s).strip())

def list_level_keys(series, it):
    out = []
    for name in series.iterations[it].meshes:
        m = REC_NAME_RE.match(name)
        if m:
            out.append((int(m.group(2)), int(m.group(1)), name))  # (lvl, patch, key)
    # process from coarse -> fine so finer overwrites later
    out.sort(key=lambda t: (t[0], t[1]))
    return out

def get_spacing_offset(mesh):
    sp  = np.array(mesh.get_attribute("gridSpacing"), dtype=float)        # (dz, dy, dx)
    off = np.array(mesh.get_attribute("gridGlobalOffset"), dtype=float)   # (z0, y0, x0)
    return sp, off

def central_index(off, sp, size, axis, plane_value=0.0, cell_centered=True):
    ax = {"z":0, "y":1, "x":2}[axis]
    shift = 0.5 if cell_centered else 0.0
    idx = int(round((plane_value - off[ax]) / sp[ax] - shift))
    return max(0, min(size[ax]-1, idx))

def read_plane(series, it, mesh_key, axis="z"):
    """Read a single-cell-thick slab through 0-plane without loading full 3D."""
    itobj = series.iterations[it]
    mesh  = itobj.meshes[mesh_key]
    rc    = mesh[REC_COMP]
    nz, ny, nx = [int(s) for s in rc.shape]
    sp, off    = get_spacing_offset(mesh)

    if axis == "z":
        iz   = central_index(off, sp, (nz,ny,nx), "z", plane_value=0.0)
        slab = rc[iz:iz+1, :, :]
        series.flush()
        arr  = np.asarray(slab)[0, :, :]
        y = off[1] + (np.arange(ny)+0.5)*sp[1]
        x = off[2] + (np.arange(nx)+0.5)*sp[2]
        y = y * CU_TO_KM
        x = x * CU_TO_KM
        return arr, (y, x)

    if axis == "y":
        iy   = central_index(off, sp, (nz,ny,nx), "y", plane_value=0.0)
        slab = rc[:, iy:iy+1, :]
        series.flush()
        arr  = np.asarray(slab)[:, 0, :]
        z = off[0] + (np.arange(nz)+0.5)*sp[0]
        x = off[2] + (np.arange(nx)+0.5)*sp[2]
        z = z * CU_TO_KM
        x = x * CU_TO_KM
        return arr, (z, x)

    raise ValueError("axis must be 'z' or 'y'")

# ---- NEW: tiny binary erosion (no SciPy needed) ----
def _erode(mask: np.ndarray, n: int = 1) -> np.ndarray:
    """Cheap binary erosion without SciPy: remove n-pixel boundary from True regions."""
    m = mask.astype(bool).copy()
    for _ in range(max(0, int(n))):
        up    = np.zeros_like(m); up[1:,  :]  = m[:-1, :]
        down  = np.zeros_like(m); down[:-1, :] = m[1:,  :]
        left  = np.zeros_like(m); left[:, 1:]  = m[:, :-1]
        right = np.zeros_like(m); right[:, :-1] = m[:, 1:]
        m = m & up & down & left & right
    return m

def composite_axis(series, it, axis="z", Ncanvas=NXNY):
    """
    Interpolate each patch in *linear code units*, composite onto a uniform canvas.
    Fine levels overwrite coarse only in their *interior*. A thin EDGE_FILL_PIX
    ring at patch borders is left to the already-present coarse data to kill seams.
    """
    if axis == "z":
        gy = np.linspace(YMIN, YMAX, Ncanvas)  # y in km
        gx = np.linspace(XMIN, XMAX, Ncanvas)  # x in km
    else:
        gy = np.linspace(ZMIN, ZMAX, Ncanvas)  # z in km
        gx = np.linspace(XMIN, XMAX, Ncanvas)  # x in km

    # start with NaNs, fill coarse → fine
    canvas = np.full((gy.size, gx.size), np.nan, dtype=np.float64)

    # process from coarse to fine (your list_level_keys already sorts like that)
    for lvl, patch, key in list_level_keys(series, it):
        try:
            slab, (a1, a2) = read_plane(series, it, key, axis=axis)
        except Exception as e:
            print(f"  Skip {key}: {e}")
            continue

        # quick reject if patch outside view
        if a2.max() < gx.min() or a2.min() > gx.max():   # x outside
            continue
        if a1.max() < gy.min() or a1.min() > gy.max():   # y/z outside
            continue

        # compute overlap on canvas
        j0 = np.searchsorted(gy, max(a1.min(), gy.min()), side="left")
        j1 = np.searchsorted(gy, min(a1.max(), gy.max()), side="right")
        i0 = np.searchsorted(gx, max(a2.min(), gx.min()), side="left")
        i1 = np.searchsorted(gx, min(a2.max(), gx.max()), side="right")
        if j1 <= j0 or i1 <= i0:
            continue

        subGy = gy[j0:j1]
        subGx = gx[i0:i1]

        # ---- interpolate in *linear* code units (as requested) ----
        vals = np.where(slab > 0.0, slab, np.nan)
        if HAVE_SCIPY:
            try:
                f = RegularGridInterpolator((a1, a2), vals,
                                            bounds_error=False, fill_value=np.nan)
                XX, YY = np.meshgrid(subGx, subGy, indexing="xy")
                sub = f((YY, XX))  # shape (subGy, subGx)
            except Exception:
                # fallback: nearest-like via index picking
                jj = np.clip(np.searchsorted(a1, subGy), 0, len(a1)-1)
                ii = np.clip(np.searchsorted(a2, subGx), 0, len(a2)-1)
                sub = vals[jj[:, None], ii[None, :]]
        else:
            jj = np.clip(np.searchsorted(a1, subGy), 0, len(a1)-1)
            ii = np.clip(np.searchsorted(a2, subGx), 0, len(a2)-1)
            sub = vals[jj[:, None], ii[None, :]]

        # validity mask for this patch
        valid = np.isfinite(sub)

        # ---- key trick: keep only the *interior* of the fine patch ----
        # anything inside EDGE_FILL_PIX of the boundary stays as-is on the canvas
        core = _erode(valid, n=EDGE_FILL_PIX)

        # write coarse → fine: coarse already there; we now overwrite only core pixels
        write_mask = core
        if np.any(write_mask):
            block = canvas[j0:j1, i0:i1]
            block[write_mask] = sub[write_mask]
            canvas[j0:j1, i0:i1] = block

        # for the outer ring (~seam zone), do nothing: existing (coarser) values remain

    return canvas, (gy, gx)

def get_time_code_units(series, it):
    itobj = series.iterations[it]
    # try standard places
    try:
        return float(getattr(itobj, "time"))
    except Exception:
        pass
    try:
        return float(itobj.get_attribute("time"))
    except Exception:
        pass
    try:
        for name in itobj.meshes:
            mesh = itobj.meshes[name]
            for getter in (mesh.get_attribute, mesh[REC_COMP].get_attribute):
                try:
                    return float(getter("time"))
                except Exception:
                    pass
    except Exception:
        pass
    return None

# ------------------ plotting ------------------
def make_grid():
    fig = plt.figure(figsize=(7.8, 3.5), constrained_layout=True)
    grid = ImageGrid(fig, 111,
                     nrows_ncols=(1, 2),
                     axes_pad=(0.45, 0.25),
                     share_all=False,
                     cbar_mode="single",
                     cbar_location="right",
                     cbar_size="3%",
                     cbar_pad=0.05,
                     label_mode="all")
    for ax in grid:
        ax.tick_params(labelleft=True)
    return fig, grid

def plot_panel(ax, data2d, extent, y_label):
    var = _clean_var(y_label)
    im = ax.imshow(
        data2d, origin="lower", extent=extent,
        norm=LogNorm(vmin=RHO_VMIN, vmax=RHO_VMAX),
        cmap=CMAP, interpolation="none"  # keep NaN boundaries crisp (no speckle)
    )
    ax.set_xlabel(r"$x$ [km]")
    ax.set_ylabel(rf"${var}$ [km]")
    minor = ticker.MultipleLocator(10)
    major = ticker.MultipleLocator(50)
    ax.xaxis.set_minor_locator(minor); ax.yaxis.set_minor_locator(minor)
    ax.xaxis.set_major_locator(major); ax.yaxis.set_major_locator(major)
    return im

def gather_series_files(path):
    return sorted(f for f in glob.glob(os.path.join(path, "*.bp*"))
                  if ".md." not in f and not f.endswith(".dir"))

def parse_itnum(fname):
    if ".it" in fname:
        tail = fname.split(".it")[-1]
        for sep in (".bp5", ".bp"):
            if sep in tail:
                try:
                    return int(tail.split(sep)[0])
                except Exception:
                    pass
    return 0

# ------------------ main ------------------
def main():
    files = gather_series_files(DATA_DIR)
    if not files:
        print("No .bp* files in", DATA_DIR); return
    frames = []

    for f in files:
        it = parse_itnum(f)
        print(f"Loading iteration {it}")
        series = io.Series(f, io.Access.read_only)

        rho_xy_cu, _ = composite_axis(series, it, axis="z", Ncanvas=NXNY)  # (y,x)
        rho_xz_cu, _ = composite_axis(series, it, axis="y", Ncanvas=NXNY)  # (z,x)

        # convert CU → cgs AFTER interpolation/compositing
        rho_xy = rho_xy_cu * RHO_CU_TO_CGS
        rho_xz = rho_xz_cu * RHO_CU_TO_CGS

        fig, grid = make_grid()

        im1 = plot_panel(grid[0], rho_xy, EXTENT_XY, y_label="y")
        im2 = plot_panel(grid[1], rho_xz, EXTENT_XZ, y_label="z")

        cbar = grid.cbar_axes[0].colorbar(im1)
        cbar.ax.set_ylabel(r"$\rho$ [g/cm$^3$]")

        t_cu = get_time_code_units(series, it)
        if t_cu is not None:
            t_ms = t_cu * TIME_CU_TO_MS
            # ---- NEW: shifted display time ----
            t_rel = t_ms - MERGER_TIME_MS

            # ---- NEW: stop when beyond desired post-merger window ----
            if (FINAL_AFTER_MS is not None) and (t_rel > FINAL_AFTER_MS):
                plt.close(fig)
                series.close()
                break

            grid[0].text(0.98, 0.98, rf"$t = {t_rel:.1f}\,\mathrm{{ms}}$", color="w",
                         transform=grid[0].transAxes, ha="right", va="top")

        frame = os.path.join(OUT_DIR, f"frame_{it:08d}.png")
        plt.savefig(frame, bbox_inches="tight", pad_inches=0.05, dpi=300)
        plt.close(fig)
        frames.append(frame)
        series.close()

    print("Combining frames…")
    imgs = [imageio.imread(p) for p in frames]

    # Try MP4 via imageio-ffmpeg; fall back to GIF if unavailable.
    have_ffmpeg = shutil.which("ffmpeg") is not None
    mp4_ok = False
    try:
        # imageio will use the FFmpeg plugin if available (from imageio-ffmpeg)
        imageio.mimsave(MOVIE_FILE, imgs, fps=FPS)
        mp4_ok = True
    except Exception as e:
        print(f"⚠️ MP4 writer not available ({e}). Falling back to GIF…")

    if not mp4_ok:
        gif_file = os.path.splitext(MOVIE_FILE)[0] + ".gif"
        imageio.mimsave(gif_file, imgs, fps=FPS)
        print(f"✅ GIF saved to: {gif_file}")
        if have_ffmpeg:
            print("Tip: convert GIF → MP4 with:")
            print(f"  ffmpeg -y -r {FPS} -i {os.path.join(OUT_DIR, 'frame_%08d.png')} "
                  f"-c:v libx264 -pix_fmt yuv420p -crf 20 {MOVIE_FILE}")
        else:
            print("Tip: install FFmpeg for direct MP4 output: pip install --user 'imageio[ffmpeg]'")

    print(f"Done. Wrote {'MP4' if mp4_ok else 'GIF'}.")

if __name__ == "__main__":
    main()

