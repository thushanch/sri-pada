"""
Sri Pada visibility study - Step 4: solar geometry.

Two products.

A) SUN-BEHIND-SUMMIT ALIGNMENT
   For every place that can see the peak, on which dates does the rising sun
   emerge from directly behind it?

   The summit sits at a fixed point in that observer's sky: bearing A, apparent
   altitude theta (both already computed by step 2). The sun passes through
   that exact point when its declination satisfies

       sin(dec) = sin(phi) sin(h) + cos(phi) cos(h) cos(A)

   with phi the observer latitude and h the sun's *true* altitude, i.e. theta
   minus astronomical refraction. Solar declination sweeps -23.44 to +23.44 and
   back once a year, so a solvable cell gets exactly two dates.

   Because sin(H) = -cos(h) sin(A) / cos(dec), any observer with A between 0 and
   180 degrees - that is, anyone west of the peak - meets the sun on its rising
   branch. Observers east of the peak never can: for them this is a sunset
   alignment instead. That single inequality is what carves the corridor out of
   the map.

B) SHADOW OF THE PEAK
   The triangular shadow Sri Pada throws at sunrise, projected onto the terrain.
   Cast by marching the DEM along the solar bearing and keeping, for each cell,
   the running maximum of z - s*tan(alt) ahead of it; a cell is shadowed if that
   maximum exceeds its own value. The index of the blocking cell is carried
   through too, so the summit's own shadow can be separated from the shadow of
   every other ridge.
"""
import os
import numpy as np
from osgeo import gdal, osr
from scipy.ndimage import map_coordinates

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
VS = os.path.join(ROOT, "Viewshed")
OUT = os.path.join(ROOT, "Solar")
DEMDIR = os.path.join(ROOT, "DEM")

SUMMIT_X, SUMMIT_Y = 469682.5, 478848.2
SUMMIT_Z = 2192.0
K_REFRAC, R_EARTH = 0.13, 6371008.8
YEAR = 2026
MAX_TIER = 150.0          # only solve where the peak is reachable at all
CO = ["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"]


# ----------------------------------------------------------------- solar core
def julian_day(y, m, d, hour=0.0):
    if m <= 2:
        y, m = y - 1, m + 12
    a = y // 100
    b = 2 - a + a // 4
    return (int(365.25 * (y + 4716)) + int(30.6001 * (m + 1))
            + d + b - 1524.5 + hour / 24.0)


def solar_declination(jd):
    """NOAA solar position, declination in degrees."""
    t = (jd - 2451545.0) / 36525.0
    l0 = (280.46646 + t * (36000.76983 + t * 0.0003032)) % 360.0
    m = 357.52911 + t * (35999.05029 - 0.0001537 * t)
    mr = np.radians(m)
    c = (np.sin(mr) * (1.914602 - t * (0.004817 + 0.000014 * t))
         + np.sin(2 * mr) * (0.019993 - 0.000101 * t)
         + np.sin(3 * mr) * 0.000289)
    true_long = l0 + c
    omega = 125.04 - 1934.136 * t
    lam = true_long - 0.00569 - 0.00478 * np.sin(np.radians(omega))
    e0 = (23.0 + (26.0 + ((21.448 - t * (46.815 + t * (0.00059 - t * 0.001813))))
                  / 60.0) / 60.0)
    e = e0 + 0.00256 * np.cos(np.radians(omega))
    return np.degrees(np.arcsin(np.sin(np.radians(e))
                                * np.sin(np.radians(lam))))


def refraction_deg(app_alt_deg):
    """Bennett: astronomical refraction from apparent altitude, degrees."""
    a = np.maximum(app_alt_deg, -0.5)
    return (1.02 / np.tan(np.radians(a + 10.3 / (a + 5.11)))) / 60.0


def declination_series(year):
    doy = np.arange(1, 366, dtype=np.float64)
    jd = np.array([julian_day(year, 1, 1, 0.5) + (d - 1) for d in doy])
    return doy, solar_declination(jd)          # ~06:00 Sri Lanka time


# ------------------------------------------------------------- lat/lon helper
def latlon_grids(gt, nx, ny, step=64):
    """Coarse exact transform, bilinear in between - sub-metre over this area."""
    sld = osr.SpatialReference(); sld.ImportFromEPSG(5235)
    sld.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    wgs = osr.SpatialReference(); wgs.ImportFromEPSG(4326)
    wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(sld, wgs)
    cs = np.arange(0, nx + step, step, dtype=np.float64)
    rs = np.arange(0, ny + step, step, dtype=np.float64)
    X = gt[0] + (cs[None, :] + 0.5) * gt[1]
    Y = gt[3] + (rs[:, None] + 0.5) * gt[5]
    X = np.repeat(X, len(rs), axis=0); Y = np.repeat(Y, len(cs), axis=1)
    pts = tr.TransformPoints(np.column_stack([X.ravel(), Y.ravel()]).tolist())
    pts = np.asarray(pts)
    lon = pts[:, 0].reshape(X.shape); lat = pts[:, 1].reshape(X.shape)
    return lat, lon, step


# =============================================================== A) alignment
def alignment(prefix, dem_path):
    print(f"\n=== alignment: {prefix} ===", flush=True)
    hr = gdal.Open(os.path.join(VS, f"{prefix}_required_height.tif"))
    an = gdal.Open(os.path.join(VS, f"{prefix}_angular_height_deg.tif"))
    dem = gdal.Open(dem_path)
    gt, wkt = hr.GetGeoTransform(), hr.GetProjection()
    nx, ny = hr.RasterXSize, hr.RasterYSize
    px = gt[1]

    latc, lonc, step = latlon_grids(gt, nx, ny)
    doy, dec = declination_series(YEAR)
    imax = int(np.argmax(dec)); imin = int(np.argmin(dec))
    # rising branch (Jan -> Jun solstice) and falling branch (Jun -> Dec)
    r_dec, r_doy = dec[:imax + 1], doy[:imax + 1]
    f_dec, f_doy = dec[imax:imin + 1][::-1], doy[imax:imin + 1][::-1]

    drv = gdal.GetDriverByName("GTiff")
    outs = {}
    for nm in ("align_doy1", "align_doy2", "align_bearing_deg"):
        o = drv.Create(os.path.join(OUT, f"{prefix}_{nm}.tif"), nx, ny, 1,
                       gdal.GDT_Float32, options=CO)
        o.SetGeoTransform(gt); o.SetProjection(wkt)
        o.GetRasterBand(1).SetNoDataValue(-9999.0)
        o.SetMetadata({"YEAR": str(YEAR),
                       "MEANING": "day-of-year the rising sun emerges from "
                                  "behind the Sri Pada summit",
                       "SUN_DISC_SEMIDIAMETER_DEG": "0.267"})
        outs[nm] = o

    blk = 512
    nsolved = 0
    bmin, bmax, dmin, dmax = 360.0, 0.0, 90.0, -90.0
    for r0 in range(0, ny, blk):
        r1 = min(r0 + blk, ny)
        h = hr.GetRasterBand(1).ReadAsArray(0, r0, nx, r1 - r0)
        th = an.GetRasterBand(1).ReadAsArray(0, r0, nx, r1 - r0)

        yy = gt[3] + (np.arange(r0, r1, dtype=np.float64)[:, None] + 0.5) * gt[5]
        xx = gt[0] + (np.arange(nx, dtype=np.float64)[None, :] + 0.5) * gt[1]
        dx = xx - SUMMIT_X
        dy = yy - SUMMIT_Y
        # arctan2(dx, dy) is summit -> observer; the sky geometry needs the
        # reciprocal, the bearing the observer looks along to find the peak
        bearing = (np.degrees(np.arctan2(dx, dy)) + 180.0) % 360.0

        rr = (np.arange(r0, r1, dtype=np.float64)[:, None]) / step
        cc = (np.arange(nx, dtype=np.float64)[None, :]) / step
        lat = map_coordinates(latc, [np.broadcast_to(rr, bearing.shape).ravel(),
                                     np.broadcast_to(cc, bearing.shape).ravel()],
                              order=1, mode="nearest").reshape(bearing.shape)

        # apparent altitude of the summit; for cells reachable only from a
        # rooftop, recompute theta at the required height rather than eye level
        d = np.hypot(dx, dy)
        cd = (1.0 - K_REFRAC) * d * d / (2.0 * R_EARTH)
        ok = (h != -9999.0) & (h <= MAX_TIER)
        theta = np.where(th != -9999.0, th, np.nan)
        need = ok & ~np.isfinite(theta)
        if need.any():
            # rooftop-only cells: put the observer at exactly the height that
            # clears the ridge, so the summit sits on their skyline
            zg = dem.GetRasterBand(1).ReadAsArray(0, r0, nx, r1 - r0)
            z_req = zg + np.maximum(h, 0.0)
            theta_roof = np.degrees(np.arctan2(SUMMIT_Z - cd - z_req,
                                               np.maximum(d, 1.0)))
            theta = np.where(need, theta_roof, theta)

        # sun must sit where the summit sits, so convert to true altitude
        h_true = theta - refraction_deg(theta)
        phi = np.radians(lat)
        A = np.radians(bearing)
        s = (np.sin(phi) * np.sin(np.radians(h_true))
             + np.cos(phi) * np.cos(np.radians(h_true)) * np.cos(A))
        # rising branch only: observer must lie west of the peak
        valid = ok & np.isfinite(theta) & (np.abs(s) <= 1.0) \
            & (bearing > 0.0) & (bearing < 180.0)
        decl = np.degrees(np.arcsin(np.clip(s, -1.0, 1.0)))
        valid &= (decl >= r_dec.min()) & (decl <= r_dec.max())

        d1 = np.where(valid, np.interp(decl, r_dec, r_doy), -9999.0)
        d2 = np.where(valid, np.interp(decl, f_dec, f_doy), -9999.0)
        nsolved += int(valid.sum())
        if valid.any():
            bv = bearing[valid]
            bmin = min(bmin, float(bv.min())); bmax = max(bmax, float(bv.max()))
            dmin = min(dmin, float(decl[valid].min()))
            dmax = max(dmax, float(decl[valid].max()))

        outs["align_doy1"].GetRasterBand(1).WriteArray(
            d1.astype(np.float32), 0, r0)
        outs["align_doy2"].GetRasterBand(1).WriteArray(
            d2.astype(np.float32), 0, r0)
        outs["align_bearing_deg"].GetRasterBand(1).WriteArray(
            np.where(ok, bearing, -9999.0).astype(np.float32), 0, r0)

    for o in outs.values():
        o.FlushCache()
    cell_km2 = (px * px) / 1e6
    print(f"  alignment solved for {nsolved:,} cells "
          f"({nsolved*cell_km2:,.0f} km2)", flush=True)
    if nsolved:
        print(f"  observer->summit bearing window {bmin:.1f} to {bmax:.1f} deg"
              f"   (theory: 66.4 to 113.6)", flush=True)
        print(f"  solar declination range {dmin:+.2f} to {dmax:+.2f} deg",
              flush=True)


# ================================================================== B) shadow
def sun_altaz(jd, lat, lon):
    """Apparent altitude/azimuth of the sun for a scalar time and location."""
    t = (jd - 2451545.0) / 36525.0
    dec = solar_declination(jd)
    l0 = (280.46646 + t * (36000.76983 + t * 0.0003032)) % 360.0
    m = 357.52911 + t * (35999.05029 - 0.0001537 * t)
    e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t)
    e0 = 23.0 + (26.0 + (21.448 - t * (46.815 + t * (0.00059 - t * 0.001813)))
                 / 60.0) / 60.0
    omega = 125.04 - 1934.136 * t
    ec = e0 + 0.00256 * np.cos(np.radians(omega))
    y = np.tan(np.radians(ec / 2.0)) ** 2
    eqt = np.degrees(
        y * np.sin(2 * np.radians(l0))
        - 2 * e * np.sin(np.radians(m))
        + 4 * e * y * np.sin(np.radians(m)) * np.cos(2 * np.radians(l0))
        - 0.5 * y * y * np.sin(4 * np.radians(l0))
        - 1.25 * e * e * np.sin(2 * np.radians(m))) * 4.0
    minutes = (jd - int(jd - 0.5) - 0.5) * 1440.0
    tst = (minutes + eqt + 4.0 * lon) % 1440.0
    ha = tst / 4.0 - 180.0
    phi, dr, hr_ = np.radians(lat), np.radians(dec), np.radians(ha)
    z = np.arccos(np.sin(phi) * np.sin(dr)
                  + np.cos(phi) * np.cos(dr) * np.cos(hr_))
    alt = 90.0 - np.degrees(z)
    # azimuth from north, eastward. The denominator is
    # sin(dec)cos(lat) - cos(dec)sin(lat)cos(H); negating it silently returns
    # 180-A, which mirrors every sunrise about due east.
    az = np.degrees(np.arctan2(-np.sin(hr_) * np.cos(dr),
                               np.sin(dr) * np.cos(phi)
                               - np.cos(hr_) * np.cos(dr) * np.sin(phi))) % 360.0
    return alt + refraction_deg(alt), az


def shadow(dem_path, date, target_alt=1.0, half_km=70):
    """Terrain shadow at the moment the sun reaches target_alt on `date`."""
    y, m, d = date
    lat0, lon0 = 6.809444, 80.499444
    # scan the local dawn window for the moment apparent altitude hits target
    jd0 = julian_day(y, m, d, 0.0)
    hrs_utc = np.linspace(-1.5, 2.5, 2000)       # 04:00-08:00 local (UTC+5:30)
    jds = jd0 + hrs_utc / 24.0
    alts = np.array([sun_altaz(j, lat0, lon0)[0] for j in jds])
    i = int(np.argmin(np.abs(alts - target_alt)))
    jd = jds[i]
    alt, az = sun_altaz(jd, lat0, lon0)
    local = (jd - jd0) * 24.0 + 5.5
    print(f"  {y}-{m:02d}-{d:02d}  local {int(local):02d}:{int(local%1*60):02d}"
          f"  sun alt {alt:.2f} deg  az {az:.2f} deg", flush=True)

    ds = gdal.Open(dem_path)
    gt = ds.GetGeoTransform(); px = gt[1]
    n = int(half_km * 1000 / px)
    c0 = int((SUMMIT_X - gt[0]) / px); r0 = int((gt[3] - SUMMIT_Y) / px)
    c_lo, r_lo = max(c0 - n, 0), max(r0 - n, 0)
    w = min(c0 + n, ds.RasterXSize) - c_lo
    hgt = min(r0 + n, ds.RasterYSize) - r_lo
    z = ds.GetRasterBand(1).ReadAsArray(c_lo, r_lo, w, hgt).astype(np.float32)
    sub_gt = (gt[0] + c_lo * px, px, 0, gt[3] - r_lo * px, 0, -px)

    # Rotate so +s points at the sun. The running maximum only ever runs along
    # s, so each perpendicular row is independent and the sweep chunks over t -
    # which keeps this inside a few hundred MB instead of a few GB.
    ux, uy = np.sin(np.radians(az)), np.cos(np.radians(az))
    diag = int(np.hypot(w, hgt)) + 2
    tt = np.arange(-diag // 2, diag // 2, dtype=np.float32)
    ss = np.arange(-diag // 2, diag // 2, dtype=np.float32)
    nS = len(ss)
    cx, cy = w / 2.0, hgt / 2.0
    sdist = ss * px
    curv = (1.0 - K_REFRAC) * sdist * sdist / (2.0 * R_EARTH)
    slope = sdist * np.tan(np.radians(alt))
    ar = np.arange(nS, dtype=np.int32)

    combo = np.zeros((len(tt), nS), np.uint8)
    CH = 256
    for t0 in range(0, len(tt), CH):
        t1 = min(t0 + CH, len(tt))
        T = tt[t0:t1][:, None]
        col = cx + T * uy + ss[None, :] * ux
        row = cy + T * ux - ss[None, :] * uy
        zr = map_coordinates(z, [row.ravel(), col.ravel()], order=1,
                             mode="constant", cval=-1e4).reshape(col.shape)
        wv = zr - slope[None, :] + curv[None, :]
        del zr, col, row

        rev = wv[:, ::-1]
        cmx = np.maximum.accumulate(rev, axis=1)
        sh = np.zeros(wv.shape, np.uint8)
        sh[:, :-1] = (cmx[:, ::-1][:, 1:] > wv[:, :-1] + 1e-3)

        # first index attaining the running max = the cell casting the shadow
        idx = np.where(rev == cmx, ar[None, :], 0)
        np.maximum.accumulate(idx, axis=1, out=idx)
        bS = ss[np.clip(nS - 1 - idx[:, ::-1], 0, nS - 1)]
        del rev, cmx, idx

        bx = sub_gt[0] + (cx + T * uy + bS * ux + 0.5) * px
        by = sub_gt[3] - (cy + T * ux - bS * uy + 0.5) * px
        near = np.hypot(bx - SUMMIT_X, by - SUMMIT_Y) < 900.0
        combo[t0:t1] = sh + 2 * (sh & near.astype(np.uint8))
        del wv, sh, bS, bx, by, near

    # rotate back to map orientation, in row blocks
    out = np.zeros((hgt, w), np.uint8)
    for r0 in range(0, hgt, 512):
        r1 = min(r0 + 512, hgt)
        rr = np.arange(r0, r1, dtype=np.float32)[:, None]
        cc = np.arange(w, dtype=np.float32)[None, :]
        dxp, dyp = cc - cx, rr - cy
        ti = (dxp * uy + dyp * ux) - tt[0]
        si = (dxp * ux - dyp * uy) - ss[0]
        out[r0:r1] = map_coordinates(
            combo, [np.broadcast_to(ti, (r1 - r0, w)).ravel(),
                    np.broadcast_to(si, (r1 - r0, w)).ravel()],
            order=0, mode="constant", cval=0).reshape(r1 - r0, w)
    del combo

    tag = f"{y}{m:02d}{d:02d}_alt{target_alt:g}"
    p = os.path.join(OUT, f"SriPada_shadow_{tag}.tif")
    drv = gdal.GetDriverByName("GTiff")
    o = drv.Create(p, w, hgt, 1, gdal.GDT_Byte, options=CO)
    o.SetGeoTransform(sub_gt)
    o.SetProjection(ds.GetProjection())
    o.GetRasterBand(1).SetNoDataValue(0)
    o.SetMetadata({"DATE": f"{y}-{m:02d}-{d:02d}", "SUN_ALT_DEG": f"{alt:.3f}",
                   "SUN_AZ_DEG": f"{az:.3f}",
                   "VALUES": "1 = terrain shadow, 3 = shadow cast by the "
                             "Sri Pada summit cone"})
    o.GetRasterBand(1).WriteArray(out)
    o.FlushCache()
    npk = int((out >= 2).sum())
    print(f"    peak shadow {npk * (px*px)/1e6:.1f} km2 -> {os.path.basename(p)}",
          flush=True)


if __name__ == "__main__":
    import sys
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    os.makedirs(OUT, exist_ok=True)
    dem30 = os.path.join(DEMDIR, "SriPada_DEM_30m.tif")
    if what in ("all", "align"):
        # island-wide 30 m grid (the 90 m tier and the 110 km near-field tier
        # are both retired; COP30's native posting is the working resolution)
        alignment("SriLanka_30m",
                  os.path.join(DEMDIR, "SriLanka_DEM_30m.tif"))
    if what in ("all", "shadow"):
        for date in [(YEAR, 3, 20), (YEAR, 6, 21), (YEAR, 12, 21),
                     (YEAR, 1, 15)]:
            for a in (0.5, 2.0):
                shadow(dem30, date, target_alt=a)
