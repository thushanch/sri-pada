"""
Sri Pada visibility study - Step 2: line-of-sight engine.

Radial sweep from the summit. For every cell on the island it answers, in one
pass, both halves of the question:

    can you see Sri Pada from here, and if not, how high would you have to be?

Method
------
Rays are cast from the summit on `naz` azimuths out to `rmax`, sampling the DEM
bilinearly every `dr` metres. Along each ray we keep the running maximum of the
vertical angle subtended at the summit,

    alpha(r) = ( z(r) - c(r) - Zs ) / r

where c(r) = (1-k) r^2 / (2R) is the combined earth-curvature and atmospheric
refraction drop (k = 0.13, the standard terrestrial coefficient; R = mean earth
radius). The running maximum *excluding the cell itself*, A*, is the grazing
angle of the highest obstruction between the summit and that range.

A cell at horizontal distance d then needs its surface to reach

    z_req = Zs + A* d + c(d)

for the summit to break its skyline. The signed quantity

    h_req = z_req - z_ground

is the master result: positive = metres of observer height needed (the tall
building question), negative = metres of headroom before the view is lost.
Visibility at eye height he is simply h_req <= he.

Because line of sight is reciprocal, sweeping outward from the summit is exact
and costs one pass instead of one viewshed per observer.

Caveats recorded in the header of every output: the DEM is a DSM (canopy is
included, so it both blocks realistically and lifts forested observers), and
COP30 truncates the summit to 2192 m against a surveyed 2243 m.
"""
import os
import time
import numpy as np
from osgeo import gdal
from scipy.ndimage import map_coordinates

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEMDIR = os.path.join(ROOT, "DEM")
OUT = os.path.join(ROOT, "Viewshed")

# ---------------------------------------------------------------- constants
SUMMIT_X, SUMMIT_Y = 469682.5, 478848.2   # EPSG:5235
SUMMIT_Z = 2192.0        # COP30 value at the summit cell (surveyed: 2243 m)
K_REFRAC = 0.13          # standard terrestrial refraction coefficient
R_EARTH = 6371008.8      # m
EYE = 1.7                # standing observer

# rooftop tiers, metres above ground
TIERS = [1.7, 10.0, 30.0, 60.0, 100.0, 150.0]

CO = ["COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"]


def _write(path, arr_gen, gt, wkt, nx, ny, dtype, nodata, block, meta=None):
    """Create a raster and fill it from a generator of (row0, array) blocks."""
    drv = gdal.GetDriverByName("GTiff")
    co = CO + (["PREDICTOR=2"] if dtype == gdal.GDT_Float32 else [])
    ds = drv.Create(path, nx, ny, 1, dtype, options=CO)
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    b = ds.GetRasterBand(1)
    if nodata is not None:
        b.SetNoDataValue(nodata)
    if meta:
        ds.SetMetadata(meta)
    for r0, a in arr_gen:
        b.WriteArray(a, 0, r0)
    b.FlushCache()
    ds = None


def run_tier(dem_path, prefix, rmax, naz, max_out_dist=None, az_chunk=1024,
             row_block=512):
    print(f"\n=== {prefix} ===")
    ds = gdal.Open(dem_path)
    gt = ds.GetGeoTransform()
    wkt = ds.GetProjection()
    nx, ny = ds.RasterXSize, ds.RasterYSize
    px = gt[1]
    dr = px
    nr = int(rmax / dr)
    print(f"DEM {nx}x{ny} @{px:g} m | rays {naz} x {nr} steps to {rmax/1000:g} km")

    z = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)

    r = ((np.arange(nr, dtype=np.float32) + 1.0) * dr)
    c_r = (1.0 - K_REFRAC) * r * r / (2.0 * R_EARTH)

    # ---------------------------------------------------- polar sweep
    Aprev = np.empty((naz, nr), np.float32)
    az = np.arange(naz, dtype=np.float64) * (2.0 * np.pi / naz)
    t0 = time.time()
    for a0 in range(0, naz, az_chunk):
        a1 = min(a0 + az_chunk, naz)
        print(f"  sweep {a0:6d}/{naz}  {time.time()-t0:6.1f}s", flush=True)
        sa = np.sin(az[a0:a1])[:, None]
        ca = np.cos(az[a0:a1])[:, None]
        X = SUMMIT_X + r[None, :] * sa
        Y = SUMMIT_Y + r[None, :] * ca
        col = (X - gt[0]) / px - 0.5
        row = (gt[3] - Y) / px - 0.5
        # outside the DEM box is open ocean -> 0 m
        zr = map_coordinates(z, [row, col], order=1, mode="constant", cval=0.0)
        alpha = (zr - c_r[None, :] - SUMMIT_Z) / r[None, :]
        np.maximum.accumulate(alpha, axis=1, out=alpha)
        Aprev[a0:a1, 0] = -1e9
        Aprev[a0:a1, 1:] = alpha[:, :-1]
        del X, Y, col, row, zr, alpha
    # wrap azimuth for interpolation
    Aprev = np.vstack([Aprev, Aprev[:1]])
    print("  polar sweep done")

    # ---------------------------------------------------- per-cell solve
    x0, y0 = gt[0], gt[3]
    naz_f = naz / (2.0 * np.pi)
    cap = max_out_dist if max_out_dist else rmax

    def blocks():
        for r0 in range(0, ny, row_block):
            r1 = min(r0 + row_block, ny)
            yy = y0 - (np.arange(r0, r1, dtype=np.float64)[:, None] + 0.5) * px
            xx = x0 + (np.arange(nx, dtype=np.float64)[None, :] + 0.5) * px
            dx = xx - SUMMIT_X
            dy = yy - SUMMIT_Y
            d = np.hypot(dx, dy)
            a = np.arctan2(dx, dy)          # bearing from north, clockwise
            np.mod(a, 2.0 * np.pi, out=a)
            ai = a * naz_f
            ri = np.clip(d / dr - 1.0, 0.0, nr - 1.0)
            astar = map_coordinates(Aprev, [ai.ravel(), ri.ravel()], order=1,
                                    mode="nearest").reshape(d.shape)
            cd = (1.0 - K_REFRAC) * d * d / (2.0 * R_EARTH)
            z_req = SUMMIT_Z + astar * d + cd
            zg = z[r0:r1, :]
            # The -1e9 sentinel that marks "nothing can block the first radial
            # step" would otherwise propagate into h as ~-1e10 on the handful
            # of cells sitting on the summit cone itself. Floor it: 3 km of
            # headroom and 3 km of headroom are the same answer.
            h = np.clip(z_req - zg, -3000.0, None).astype(np.float32)
            yield r0, d.astype(np.float32), h

    meta = {
        "SUMMIT_X": str(SUMMIT_X), "SUMMIT_Y": str(SUMMIT_Y),
        "SUMMIT_Z_M": str(SUMMIT_Z),
        "SUMMIT_Z_NOTE": "COP30 summit cell; surveyed height is 2243 m",
        "REFRACTION_K": str(K_REFRAC), "EARTH_R_M": str(R_EARTH),
        "DEM": "Copernicus GLO-30 DSM (canopy included)",
        "CRS": "EPSG:5235 SLD99 / Sri Lanka Grid 1999",
    }

    # We need three passes' worth of products from one computation, so cache
    # each block's result and write to all outputs as we go.
    drv = gdal.GetDriverByName("GTiff")
    outs = {}
    for name, dt, nd in (("required_height", gdal.GDT_Float32, -9999.0),
                         ("distance_km", gdal.GDT_Float32, -9999.0),
                         ("class", gdal.GDT_Byte, 255),
                         ("angular_height_deg", gdal.GDT_Float32, -9999.0)):
        p = os.path.join(OUT, f"{prefix}_{name}.tif")
        o = drv.Create(p, nx, ny, 1, dt, options=CO)
        o.SetGeoTransform(gt); o.SetProjection(wkt)
        o.GetRasterBand(1).SetNoDataValue(nd)
        o.SetMetadata(meta)
        outs[name] = o

    tier_arr = np.array(TIERS, dtype=np.float32)
    counts = np.zeros(len(TIERS) + 1, dtype=np.int64)
    ncells_in = 0

    t1 = time.time()
    for r0, d, h in blocks():
        if (r0 // row_block) % 4 == 0:
            print(f"  solve row {r0:6d}/{ny}  {time.time()-t1:6.1f}s", flush=True)
        inrange = d <= cap
        ncells_in += int(inrange.sum())

        hh = np.where(inrange, h, -9999.0).astype(np.float32)
        outs["required_height"].GetRasterBand(1).WriteArray(hh, 0, r0)

        dk = np.where(inrange, d / 1000.0, -9999.0).astype(np.float32)
        outs["distance_km"].GetRasterBand(1).WriteArray(dk, 0, r0)

        # class: index of the first tier that clears it, 250 = never
        cls = np.full(d.shape, 250, dtype=np.uint8)
        for i in range(len(TIERS) - 1, -1, -1):
            cls = np.where(h <= tier_arr[i], np.uint8(i), cls)
        cls = np.where(inrange, cls, np.uint8(255)).astype(np.uint8)
        outs["class"].GetRasterBand(1).WriteArray(cls, 0, r0)
        for i in range(len(TIERS)):
            counts[i] += int((cls == i).sum())
        counts[-1] += int((cls == 250).sum())

        # apparent elevation angle of the summit for a standing observer,
        # only where actually visible from the ground
        vis = (h <= EYE) & inrange
        zg = z[r0:r0 + d.shape[0], :]
        cd = (1.0 - K_REFRAC) * d * d / (2.0 * R_EARTH)
        with np.errstate(invalid="ignore", divide="ignore"):
            ang = np.degrees(np.arctan2(SUMMIT_Z - cd - (zg + EYE),
                                        np.maximum(d, 1.0)))
        ang = np.where(vis, ang, -9999.0).astype(np.float32)
        outs["angular_height_deg"].GetRasterBand(1).WriteArray(ang, 0, r0)

    for o in outs.values():
        o.FlushCache()
    outs = None
    del z, Aprev

    cell_km2 = (px * px) / 1e6
    print(f"  cells analysed: {ncells_in:,}  ({ncells_in*cell_km2:,.0f} km2)")
    lbl = [f"<= {t:g} m" for t in TIERS] + ["never"]
    cum = 0
    for i, n in enumerate(counts):
        if i < len(TIERS):
            cum += n
        print(f"    height needed {lbl[i]:>9s}: {n*cell_km2:10,.0f} km2"
              + (f"   cumulative {cum*cell_km2:,.0f} km2" if i < len(TIERS) else ""))
    return counts, cell_km2


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)

    # near field: 30 m, results kept inside 110 km
    run_tier(os.path.join(DEMDIR, "SriPada_DEM_30m.tif"), "SriPada_30m",
             rmax=110_000, naz=21600, max_out_dist=110_000)

    # island wide: 90 m, rays long enough to reach every corner of the grid
    run_tier(os.path.join(DEMDIR, "SriLanka_DEM_90m.tif"), "SriLanka_90m",
             rmax=400_000, naz=14400)
