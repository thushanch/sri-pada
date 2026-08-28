"""
Sri Pada visibility study - Step 11: line-of-sight engine, 30 m island-wide.

Same physics as before, rebuilt to run over 150 M cells on a modest machine.

PHYSICS (unchanged, stated in full because it is the point of the study)

  Earth curvature and atmospheric refraction are both modelled, everywhere,
  as the standard combined drop

      c(d) = (1 - k) d^2 / (2R),    k = 0.13,  R = 6,371,008.8 m

  k is the terrestrial refraction coefficient: light bends toward the earth,
  so the effective curvature is reduced to (1-k) of geometric. This is applied
  twice - once to every terrain sample along a ray, so blocking ridges sink by
  the right amount, and once at the target cell. It is not a cosmetic
  correction: it moves the sea-level horizon for a 2192 m peak from 172 km
  (k=0) to 184 km (k=0.13), and drops the summit 351 m at Colombo's range.

  Ray sweep from the summit, running maximum of the vertical angle

      alpha(r) = ( z(r) - c(r) - Zs ) / r

  gives A*, the grazing angle of the highest obstruction inside range d. The
  cell then needs its surface to reach z_req = Zs + A* d + c(d), and the signed
  result h_req = z_req - z_ground is the master layer: positive = observer
  height required, negative = headroom in hand.

MEMORY / SPEED

  The DEM is held as int16 decimetres (0.1 m quantisation, well inside COP30's
  ~2-4 m vertical accuracy) - 300 MB instead of 600 MB. The solve runs tile by
  tile rather than row by row: a tile away from the summit subtends a narrow
  azimuth wedge, so only a thin contiguous slice of the polar array is touched
  and it stays in cache. Tiles that straddle the north seam fall back to the
  full array.

  Rays stop at 190 km. Beyond that the curvature drop alone (2465 m at 190 km)
  already exceeds the summit's 2192 m, so no sea-level observer can see it at
  any range past ~184 km; those cells are written as never-visible.
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

SUMMIT_X, SUMMIT_Y = 469682.5, 478848.2
SUMMIT_Z = 2192.0
K_REFRAC = 0.13
R_EARTH = 6371008.8
EYE = 1.7

RMAX = 190_000.0
NAZ = 24000
TILE = 1024

TIERS = [1.7, 10.0, 30.0, 60.0, 100.0, 150.0]
CO = ["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"]


def main():
    t0 = time.time()
    dem_path = os.path.join(DEMDIR, "SriLanka_DEM_30m.tif")
    ds = gdal.Open(dem_path)
    gt = ds.GetGeoTransform()
    wkt = ds.GetProjection()
    nx, ny = ds.RasterXSize, ds.RasterYSize
    px = gt[1]
    dr = px
    nr = int(RMAX / dr)
    print(f"DEM {nx}x{ny} @{px:g} m | {NAZ} rays x {nr} steps to {RMAX/1000:g} km",
          flush=True)

    # ---- DEM as int16 decimetres
    print("loading DEM ...", flush=True)
    z16 = np.empty((ny, nx), np.int16)
    BL = 2048
    for r0 in range(0, ny, BL):
        h = min(BL, ny - r0)
        blk = ds.GetRasterBand(1).ReadAsArray(0, r0, nx, h)
        np.clip(blk * 10.0, -32000, 32000, out=blk)
        z16[r0:r0 + h] = blk.astype(np.int16)
    del blk
    print(f"  {z16.nbytes/1e6:.0f} MB  ({time.time()-t0:.0f}s)", flush=True)

    # ---- polar sweep
    r = ((np.arange(nr, dtype=np.float32) + 1.0) * dr)
    c_r = (1.0 - K_REFRAC) * r * r / (2.0 * R_EARTH)
    Aprev = np.empty((NAZ + 1, nr), np.float32)
    az = np.arange(NAZ, dtype=np.float64) * (2.0 * np.pi / NAZ)
    CH = 512
    t1 = time.time()
    for a0 in range(0, NAZ, CH):
        a1 = min(a0 + CH, NAZ)
        sa = np.sin(az[a0:a1])[:, None]
        ca = np.cos(az[a0:a1])[:, None]
        col = (SUMMIT_X + r[None, :] * sa - gt[0]) / px - 0.5
        row = (gt[3] - (SUMMIT_Y + r[None, :] * ca)) / px - 0.5
        zr = map_coordinates(z16, [row, col], order=1, mode="constant",
                             cval=0.0, output=np.float32)
        zr *= 0.1                                    # decimetres -> metres
        alpha = (zr - c_r[None, :] - SUMMIT_Z) / r[None, :]
        np.maximum.accumulate(alpha, axis=1, out=alpha)
        Aprev[a0:a1, 0] = -10.0        # nothing can block the first step
        Aprev[a0:a1, 1:] = alpha[:, :-1]
        del col, row, zr, alpha
        if (a0 // CH) % 8 == 0:
            print(f"  sweep {a0:6d}/{NAZ}  {time.time()-t1:6.1f}s", flush=True)
    Aprev[NAZ] = Aprev[0]              # wrap
    print(f"  polar sweep done, {Aprev.nbytes/1e6:.0f} MB "
          f"({time.time()-t1:.0f}s)", flush=True)

    # ---- outputs
    drv = gdal.GetDriverByName("GTiff")
    outs = {}
    for name, dt, nd in (("required_height", gdal.GDT_Float32, -9999.0),
                         ("class", gdal.GDT_Byte, 255),
                         ("angular_height_deg", gdal.GDT_Float32, -9999.0),
                         ("distance_km", gdal.GDT_Float32, -9999.0)):
        p = os.path.join(OUT, f"SriLanka_30m_{name}.tif")
        o = drv.Create(p, nx, ny, 1, dt, options=CO)
        o.SetGeoTransform(gt)
        o.SetProjection(wkt)
        o.GetRasterBand(1).SetNoDataValue(nd)
        o.SetMetadata({
            "SUMMIT_XY": f"{SUMMIT_X},{SUMMIT_Y}", "SUMMIT_Z_M": str(SUMMIT_Z),
            "SUMMIT_Z_NOTE": "COP30 summit cell; surveyed 2243 m",
            "REFRACTION_K": str(K_REFRAC), "EARTH_R_M": str(R_EARTH),
            "CURV_REFRAC_MODEL": "c(d)=(1-k)d^2/(2R), applied to every ray "
                                 "sample and to the target cell",
            "RESOLUTION_M": "30", "RMAX_M": str(int(RMAX)),
            "DEM": "Copernicus GLO-30 DSM (canopy included)",
            "CRS": "EPSG:5235 SLD99 / Sri Lanka Grid 1999",
        })
        outs[name] = o

    tier = np.array(TIERS, np.float32)
    naz_f = NAZ / (2.0 * np.pi)
    counts = np.zeros(len(TIERS) + 1, np.int64)
    t2 = time.time()
    ntile = ((ny + TILE - 1) // TILE) * ((nx + TILE - 1) // TILE)
    done = 0

    for ry in range(0, ny, TILE):
        th = min(TILE, ny - ry)
        yy = gt[3] - (np.arange(ry, ry + th, dtype=np.float64)[:, None] + 0.5) * px
        for rx in range(0, nx, TILE):
            tw = min(TILE, nx - rx)
            xx = gt[0] + (np.arange(rx, rx + tw, dtype=np.float64)[None, :] + 0.5) * px
            dx = xx - SUMMIT_X
            dy = yy - SUMMIT_Y
            d = np.hypot(dx, dy)
            a = np.arctan2(dx, dy)
            np.mod(a, 2.0 * np.pi, out=a)
            ai = a * naz_f
            ri = np.clip(d / dr - 1.0, 0.0, nr - 1.0)

            # narrow azimuth wedge -> touch only a thin slice of the sweep
            lo, hi = ai.min(), ai.max()
            if hi - lo > NAZ * 0.5:
                sub, off = Aprev, 0.0
            else:
                i0 = max(int(np.floor(lo)) - 1, 0)
                i1 = min(int(np.ceil(hi)) + 2, NAZ + 1)
                sub, off = Aprev[i0:i1], float(i0)

            astar = map_coordinates(sub, [(ai - off).ravel(), ri.ravel()],
                                    order=1, mode="nearest").reshape(d.shape)
            cd = (1.0 - K_REFRAC) * d * d / (2.0 * R_EARTH)
            zg = z16[ry:ry + th, rx:rx + tw].astype(np.float32) * 0.1
            h = np.clip(SUMMIT_Z + astar * d + cd - zg, -3000.0, None)

            inr = d <= RMAX
            hh = np.where(inr, h, 30000.0).astype(np.float32)
            outs["required_height"].GetRasterBand(1).WriteArray(hh, rx, ry)
            outs["distance_km"].GetRasterBand(1).WriteArray(
                (d / 1000.0).astype(np.float32), rx, ry)

            cls = np.full(d.shape, 250, np.uint8)
            for i in range(len(TIERS) - 1, -1, -1):
                cls = np.where(hh <= tier[i], np.uint8(i), cls)
            outs["class"].GetRasterBand(1).WriteArray(cls, rx, ry)
            for i in range(len(TIERS)):
                counts[i] += int((cls == i).sum())
            counts[-1] += int((cls == 250).sum())

            vis = hh <= EYE
            ang = np.degrees(np.arctan2(SUMMIT_Z - cd - (zg + EYE),
                                        np.maximum(d, 1.0)))
            outs["angular_height_deg"].GetRasterBand(1).WriteArray(
                np.where(vis, ang, -9999.0).astype(np.float32), rx, ry)

            done += 1
            if done % 40 == 0:
                el = time.time() - t2
                print(f"  tile {done}/{ntile}  {el:6.0f}s  "
                      f"eta {el/done*(ntile-done):6.0f}s", flush=True)

    for o in outs.values():
        o.FlushCache()
    outs = None

    cell = (px * px) / 1e6
    print(f"\nsolve done in {time.time()-t2:.0f}s "
          f"(total {time.time()-t0:.0f}s)", flush=True)
    lbl = [f"<= {t:g} m" for t in TIERS] + ["never"]
    cum = 0
    for i, n in enumerate(counts):
        if i < len(TIERS):
            cum += n
        print(f"  {lbl[i]:>9s}: {n*cell:10,.0f} km2"
              + (f"   cum {cum*cell:,.0f}" if i < len(TIERS) else ""))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    main()
