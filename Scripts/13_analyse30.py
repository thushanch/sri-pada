"""
Sri Pada visibility study - Step 13: statistics and point results, 30 m.

Everything block-wise: the 30 m rasters are 150 M cells and will not fit in
RAM alongside each other on this machine.
"""
import csv
import datetime as dt
import json
import os
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEMDIR, VS, OSMD, SOL = (os.path.join(ROOT, d)
                         for d in ("DEM", "Viewshed", "OSM", "Solar"))
OUTD = os.path.join(ROOT, "Results")
os.makedirs(OUTD, exist_ok=True)

SUMMIT_X, SUMMIT_Y, SUMMIT_Z = 469682.5, 478848.2, 2192.0
K, RE, EYE = 0.13, 6371008.8, 1.7
TIERS = [1.7, 10.0, 30.0, 60.0, 100.0, 150.0]
TIER_LBL = ["ground level", "<=10 m (3 storeys)", "<=30 m (10 storeys)",
            "<=60 m (20 storeys)", "<=100 m (33 storeys)",
            "<=150 m (50 storeys)"]
GRADES = [(2.0, 1, "prominent"), (1.0, 2, "clear"), (0.5, 3, "distinct"),
          (0.2, 4, "faint"), (0.0, 5, "marginal")]
CO = ["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES"]
BL = 1024

P_CLS = os.path.join(VS, "SriLanka_30m_class.tif")
P_HR = os.path.join(VS, "SriLanka_30m_required_height.tif")
P_ANG = os.path.join(VS, "SriLanka_30m_angular_height_deg.tif")
P_LM = os.path.join(DEMDIR, "SriLanka_landmask_30m.tif")


def stats_and_grade():
    cls_ds, lm_ds, ang_ds = gdal.Open(P_CLS), gdal.Open(P_LM), gdal.Open(P_ANG)
    nx, ny = cls_ds.RasterXSize, cls_ds.RasterYSize
    gt, wkt = cls_ds.GetGeoTransform(), cls_ds.GetProjection()
    cell = (gt[1] * abs(gt[5])) / 1e6

    g_out = gdal.GetDriverByName("GTiff").Create(
        os.path.join(VS, "SriLanka_30m_grade.tif"), nx, ny, 1,
        gdal.GDT_Byte, options=CO)
    g_out.SetGeoTransform(gt); g_out.SetProjection(wkt)
    g_out.GetRasterBand(1).SetNoDataValue(0)
    g_out.SetMetadata({"1": "prominent >=2 deg", "2": "clear 1-2",
                       "3": "distinct 0.5-1", "4": "faint 0.2-0.5",
                       "5": "marginal <0.2"})

    counts = np.zeros(7, np.int64)
    gcounts = np.zeros(6, np.int64)
    land_total = 0
    for r0 in range(0, ny, BL):
        h = min(BL, ny - r0)
        c = cls_ds.GetRasterBand(1).ReadAsArray(0, r0, nx, h)
        l = lm_ds.GetRasterBand(1).ReadAsArray(0, r0, nx, h)
        a = ang_ds.GetRasterBand(1).ReadAsArray(0, r0, nx, h)
        m = l == 1
        land_total += int(m.sum())
        for i in range(6):
            counts[i] += int(((c == i) & m).sum())
        counts[6] += int(((c == 250) & m).sum())
        g = np.zeros(c.shape, np.uint8)
        for thr, code, _ in GRADES:
            g = np.where((a >= thr) & (a != -9999.0) & (g == 0),
                         np.uint8(code), g)
        g = np.where(m, g, np.uint8(0)).astype(np.uint8)
        for _, code, _ in GRADES:
            gcounts[code] += int((g == code).sum())
        g_out.GetRasterBand(1).WriteArray(g, 0, r0)
    g_out.FlushCache(); g_out = None

    tot = land_total * cell
    print(f"\n=== 30 m island-wide, LAND ONLY ({tot:,.0f} km2) ===")
    cum = 0.0
    for i, lbl in enumerate(TIER_LBL):
        a_ = counts[i] * cell; cum += a_
        print(f"  {lbl:22s} {a_:9,.0f} km2   cum {cum:9,.0f}  ({100*cum/tot:5.1f}%)")
    print(f"  {'never visible':22s} {counts[6]*cell:9,.0f} km2"
          f"  ({100*counts[6]*cell/tot:5.1f}%)")
    print("\n  how strongly it reads (ground-level viewers):")
    for thr, code, name in GRADES:
        print(f"    {name:10s} {gcounts[code]*cell:9,.0f} km2")
    return tot


class BlockSampler:
    """Sample many scattered points without loading 150 M cells."""

    def __init__(self, paths):
        self.ds = {k: gdal.Open(p) for k, p in paths.items() if os.path.exists(p)}
        any_ = next(iter(self.ds.values()))
        self.gt = any_.GetGeoTransform()
        self.nx, self.ny = any_.RasterXSize, any_.RasterYSize

    def sample(self, xs, ys):
        n = len(xs)
        cols = ((np.asarray(xs) - self.gt[0]) / self.gt[1]).astype(np.int64)
        rows = ((self.gt[3] - np.asarray(ys)) / abs(self.gt[5])).astype(np.int64)
        ok = (cols >= 0) & (cols < self.nx) & (rows >= 0) & (rows < self.ny)
        out = {k: np.full(n, np.nan) for k in self.ds}
        order = np.argsort(rows)
        for r0 in range(0, self.ny, BL):
            sel = order[(rows[order] >= r0) & (rows[order] < r0 + BL) & ok[order]]
            if not len(sel):
                continue
            h = min(BL, self.ny - r0)
            for k, d in self.ds.items():
                arr = d.GetRasterBand(1).ReadAsArray(0, r0, self.nx, h)
                out[k][sel] = arr[rows[sel] - r0, cols[sel]]
        return out, ok


def classify(h):
    for i, t in enumerate(TIERS):
        if h <= t:
            return i, TIER_LBL[i]
    return 250, "never"


def doy_to_date(d):
    if d is None or not np.isfinite(d) or d <= 0:
        return ""
    return (dt.date(2026, 1, 1) + dt.timedelta(days=float(d) - 1)).strftime("%d %b")


if __name__ == "__main__":
    stats_and_grade()

    wgs = osr.SpatialReference(); wgs.ImportFromEPSG(4326)
    wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    sld = osr.SpatialReference(); sld.ImportFromEPSG(5235)
    sld.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(wgs, sld)

    paths = {"h": P_HR, "ang": P_ANG, "lm": P_LM,
             "d1": os.path.join(SOL, "SriLanka_30m_align_doy1.tif"),
             "d2": os.path.join(SOL, "SriLanka_30m_align_doy2.tif")}
    smp = BlockSampler(paths)

    def load(p):
        with open(p, encoding="utf-8") as f:
            return json.load(f)["features"]

    # ---------------- places
    feats = load(os.path.join(OSMD, "SriLanka_places.geojson"))
    lons = [f["geometry"]["coordinates"][0] for f in feats]
    lats = [f["geometry"]["coordinates"][1] for f in feats]
    pts = tr.TransformPoints(list(zip(lons, lats)))
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    vals, ok = smp.sample(xs, ys)

    rows = []
    for i, f in enumerate(feats):
        if not ok[i] or not np.isfinite(vals["h"][i]):
            continue
        h = float(vals["h"][i])
        code, lbl = classify(h)
        d = float(np.hypot(xs[i] - SUMMIT_X, ys[i] - SUMMIT_Y))
        brg = float(np.degrees(np.arctan2(SUMMIT_X - xs[i],
                                          SUMMIT_Y - ys[i])) % 360)
        ang = float(vals["ang"][i])
        p = f["properties"]
        rows.append({
            "name": p.get("name", ""), "place": p.get("place", ""),
            "lon": round(lons[i], 5), "lat": round(lats[i], 5),
            "distance_km": round(d / 1000, 1),
            "bearing_to_summit_deg": round(brg, 1),
            "required_height_m": round(h, 1), "visibility": lbl,
            "angular_height_deg": round(ang, 3) if ang != -9999.0 else "",
            "sunrise_align_1": doy_to_date(vals["d1"][i]) if "d1" in vals else "",
            "sunrise_align_2": doy_to_date(vals["d2"][i]) if "d2" in vals else "",
        })
    rows.sort(key=lambda r: r["distance_km"])
    fp = os.path.join(OUTD, "SriPada_places_visibility.csv")
    with open(fp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    vis = [r for r in rows if r["visibility"] == "ground level"]
    roof = [r for r in rows if r["visibility"] not in ("ground level", "never")]
    print(f"\nplaces {len(rows)}: ground {len(vis)}, rooftop-only {len(roof)}, "
          f"never {len(rows)-len(vis)-len(roof)}")

    with open(os.path.join(OUTD, "SriPada_places_visibility.geojson"), "w",
              encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": [
            {"type": "Feature", "geometry": {"type": "Point",
             "coordinates": [r["lon"], r["lat"]]}, "properties": r}
            for r in rows]}, f, ensure_ascii=False)

    # ---------------- buildings
    bf = load(os.path.join(OSMD, "SriLanka_buildings.geojson"))
    bf = [f for f in bf if f["properties"].get("height_m") is not None]
    blons = [f["geometry"]["coordinates"][0] for f in bf]
    blats = [f["geometry"]["coordinates"][1] for f in bf]
    bpts = tr.TransformPoints(list(zip(blons, blats)))
    bxs = [p[0] for p in bpts]; bys = [p[1] for p in bpts]
    bvals, bok = smp.sample(bxs, bys)

    brows = []
    for i, f in enumerate(bf):
        if not bok[i] or not np.isfinite(bvals["h"][i]):
            continue
        need = float(bvals["h"][i])
        hb = float(f["properties"]["height_m"])
        if need <= EYE or need > 300 or hb < need:
            continue
        d = float(np.hypot(bxs[i] - SUMMIT_X, bys[i] - SUMMIT_Y))
        brows.append({"name": f["properties"].get("name", ""),
                      "lon": round(blons[i], 5), "lat": round(blats[i], 5),
                      "building_height_m": hb,
                      "required_height_m": round(need, 1),
                      "spare_m": round(hb - need, 1),
                      "distance_km": round(d / 1000, 1)})
    brows.sort(key=lambda r: -r["spare_m"])
    fb = os.path.join(OUTD, "SriPada_buildings_visibility.csv")
    with open(fb, "w", newline="", encoding="utf-8-sig") as f:
        if brows:
            w = csv.DictWriter(f, fieldnames=list(brows[0].keys()))
            w.writeheader(); w.writerows(brows)
    print(f"buildings {len(bf):,} tested, roof wins the view: {len(brows):,}")
    for r in brows[:8]:
        print(f"   {(r['name'] or '(unnamed)')[:32]:34s}{r['building_height_m']:6.0f} m "
              f"needs {r['required_height_m']:6.1f} m @{r['distance_km']:6.1f} km")
