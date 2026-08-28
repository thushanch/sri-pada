"""
Sri Pada visibility study - Step 5: land masking, statistics, point results.

The raw viewshed covers everything inside the analysis box, most of which is
the Indian Ocean. Copernicus codes sea as exactly 0.0 m, which recovers the
coastline to within 1.2% of Sri Lanka's surveyed 65,610 km2, so that is what
the land mask is built from (nearest-neighbour sampled, never averaged, so the
zeros survive).

Produces:
  *_landmask.tif                land = 1
  *_grade.tif                   how strongly the peak reads, 1-5, land only
  SriPada_places_visibility.*   every OSM city/town/village, with distance,
                                bearing, required observer height, apparent
                                angular height and the two sunrise-alignment
                                dates
  SriPada_buildings_visibility.* tagged buildings whose roof clears the terrain
                                when the ground does not
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
EYE = 1.7
TIERS = [1.7, 10.0, 30.0, 60.0, 100.0, 150.0]
TIER_LBL = ["ground level", "<=10 m (3 storeys)", "<=30 m (10 storeys)",
            "<=60 m (20 storeys)", "<=100 m (33 storeys)",
            "<=150 m (50 storeys)"]
CO = ["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"]

GRADES = [(2.0, 1, "prominent"), (1.0, 2, "clear"), (0.5, 3, "distinct"),
          (0.2, 4, "faint"), (0.0, 5, "marginal")]


def build_landmask(prefix, ref_path):
    ref = gdal.Open(ref_path)
    gt, wkt = ref.GetGeoTransform(), ref.GetProjection()
    nx, ny = ref.RasterXSize, ref.RasterYSize
    bounds = (gt[0], gt[3] + ny * gt[5], gt[0] + nx * gt[1], gt[3])
    vrt = os.path.join(DEMDIR, "SriLanka_COP30_wgs84.vrt")
    tmp = gdal.Warp("", vrt, format="MEM", dstSRS="EPSG:5235",
                    xRes=gt[1], yRes=abs(gt[5]), outputBounds=bounds,
                    resampleAlg="near", outputType=gdal.GDT_Float32)
    z = tmp.GetRasterBand(1).ReadAsArray()
    land = (z != 0.0).astype(np.uint8)
    p = os.path.join(VS, f"{prefix}_landmask.tif")
    o = gdal.GetDriverByName("GTiff").Create(p, nx, ny, 1, gdal.GDT_Byte,
                                             options=CO)
    o.SetGeoTransform(gt); o.SetProjection(wkt)
    o.GetRasterBand(1).WriteArray(land)
    o.FlushCache(); o = None
    km2 = land.sum() * (gt[1] * abs(gt[5])) / 1e6
    print(f"{prefix}: land {km2:,.0f} km2")
    return p, km2


def land_stats(prefix, cap_km=None):
    cls = gdal.Open(os.path.join(VS, f"{prefix}_class.tif")).ReadAsArray()
    land = gdal.Open(os.path.join(VS, f"{prefix}_landmask.tif")).ReadAsArray()
    ds = gdal.Open(os.path.join(VS, f"{prefix}_class.tif"))
    gt = ds.GetGeoTransform()
    cell = (gt[1] * abs(gt[5])) / 1e6
    m = land == 1
    if cap_km:
        d = gdal.Open(os.path.join(VS, f"{prefix}_distance_km.tif")).ReadAsArray()
        m &= (d >= 0) & (d <= cap_km)
    tot = m.sum() * cell
    print(f"\n--- {prefix} : land only ({tot:,.0f} km2) ---")
    cum = 0.0
    rows = []
    for i, lbl in enumerate(TIER_LBL):
        a = ((cls == i) & m).sum() * cell
        cum += a
        rows.append((lbl, a, cum, 100 * cum / tot))
        print(f"  {lbl:22s} {a:9,.0f} km2   cumulative {cum:9,.0f} km2"
              f"  ({100*cum/tot:5.1f}% of land)")
    never = ((cls == 250) & m).sum() * cell
    print(f"  {'never visible':22s} {never:9,.0f} km2"
          f"  ({100*never/tot:5.1f}% of land)")
    return rows, never, tot


def build_grade(prefix):
    ang = gdal.Open(os.path.join(VS, f"{prefix}_angular_height_deg.tif"))
    a = ang.ReadAsArray()
    land = gdal.Open(os.path.join(VS, f"{prefix}_landmask.tif")).ReadAsArray()
    g = np.zeros(a.shape, np.uint8)
    for thr, code, _ in GRADES:
        g = np.where((a >= thr) & (a != -9999.0) & (g == 0), np.uint8(code), g)
    g = np.where(land == 1, g, np.uint8(0)).astype(np.uint8)
    gt, wkt = ang.GetGeoTransform(), ang.GetProjection()
    p = os.path.join(VS, f"{prefix}_grade.tif")
    o = gdal.GetDriverByName("GTiff").Create(p, ang.RasterXSize,
                                             ang.RasterYSize, 1,
                                             gdal.GDT_Byte, options=CO)
    o.SetGeoTransform(gt); o.SetProjection(wkt)
    o.GetRasterBand(1).SetNoDataValue(0)
    o.SetMetadata({"1": "prominent >=2 deg", "2": "clear 1-2 deg",
                   "3": "distinct 0.5-1 deg", "4": "faint 0.2-0.5 deg",
                   "5": "marginal <0.2 deg, needs exceptional air"})
    o.GetRasterBand(1).WriteArray(g)
    o.FlushCache(); o = None
    cell = (gt[1] * abs(gt[5])) / 1e6
    print(f"\n  visibility grade, land seeing the peak from the ground:")
    for thr, code, name in GRADES:
        print(f"    {name:10s} {(g == code).sum()*cell:9,.0f} km2")
    return p


# --------------------------------------------------------------- point sampling
class Sampler:
    def __init__(self, prefix):
        self.r = {}
        for nm in ("required_height", "distance_km", "angular_height_deg",
                   "class", "landmask"):
            p = os.path.join(VS, f"{prefix}_{nm}.tif")
            if os.path.exists(p):
                self.r[nm] = gdal.Open(p)
        for nm in ("align_doy1", "align_doy2"):
            p = os.path.join(SOL, f"{prefix}_{nm}.tif")
            if os.path.exists(p):
                self.r[nm] = gdal.Open(p)
        any_ds = next(iter(self.r.values()))
        self.gt = any_ds.GetGeoTransform()
        self.nx, self.ny = any_ds.RasterXSize, any_ds.RasterYSize
        self.arrays = {k: v.GetRasterBand(1).ReadAsArray()
                       for k, v in self.r.items()}

    def sample(self, x, y):
        c = int((x - self.gt[0]) / self.gt[1])
        r = int((self.gt[3] - y) / abs(self.gt[5]))
        if not (0 <= c < self.nx and 0 <= r < self.ny):
            return None
        return {k: float(v[r, c]) for k, v in self.arrays.items()}


def doy_to_date(doy, year=2026):
    if doy is None or doy <= 0:
        return ""
    d = dt.date(year, 1, 1) + dt.timedelta(days=float(doy) - 1)
    return d.strftime("%d %b")


def load_geojson(p):
    with open(p, encoding="utf-8") as f:
        return json.load(f)["features"]


def to_sld():
    w = osr.SpatialReference(); w.ImportFromEPSG(4326)
    w.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    s = osr.SpatialReference(); s.ImportFromEPSG(5235)
    s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return osr.CoordinateTransformation(w, s)


def classify(h):
    for i, t in enumerate(TIERS):
        if h <= t:
            return i, TIER_LBL[i]
    return 250, "never"


if __name__ == "__main__":
    print("=== land masks ===")
    build_landmask("SriPada_30m", os.path.join(VS, "SriPada_30m_class.tif"))
    build_landmask("SriLanka_90m", os.path.join(VS, "SriLanka_90m_class.tif"))

    land_stats("SriLanka_90m")
    land_stats("SriPada_30m", cap_km=110)
    build_grade("SriLanka_90m")
    build_grade("SriPada_30m")

    print("\n=== point results ===")
    s30 = Sampler("SriPada_30m")
    s90 = Sampler("SriLanka_90m")
    tr = to_sld()

    def sample_best(lon, lat):
        x, y, _ = tr.TransformPoint(lon, lat)
        v = s30.sample(x, y)
        tier = "30 m"
        if v is None or v.get("distance_km", -1) < 0:
            v = s90.sample(x, y)
            tier = "90 m"
        if v is None:
            return None, None, None, None
        return v, tier, x, y

    # ---- places
    feats = load_geojson(os.path.join(OSMD, "SriLanka_places.geojson"))
    rows = []
    for f in feats:
        lon, lat = f["geometry"]["coordinates"]
        v, tier, x, y = sample_best(lon, lat)
        if v is None:
            continue
        h = v["required_height"]
        code, lbl = classify(h)
        bearing = (np.degrees(np.arctan2(SUMMIT_X - x, SUMMIT_Y - y))) % 360
        p = f["properties"]
        rows.append({
            "name": p.get("name", ""), "name_en": p.get("name_en", ""),
            "place": p.get("place", ""), "lon": round(lon, 5),
            "lat": round(lat, 5),
            "distance_km": round(v["distance_km"], 1),
            "bearing_to_summit_deg": round(bearing, 1),
            "required_height_m": round(h, 1),
            "visibility": lbl,
            "angular_height_deg": (round(v["angular_height_deg"], 3)
                                   if v["angular_height_deg"] != -9999 else ""),
            "sunrise_align_1": doy_to_date(v.get("align_doy1")),
            "sunrise_align_2": doy_to_date(v.get("align_doy2")),
            "grid": tier,
        })
    rows.sort(key=lambda r: r["distance_km"])
    fp = os.path.join(OUTD, "SriPada_places_visibility.csv")
    with open(fp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    vis = [r for r in rows if r["visibility"] == "ground level"]
    roof = [r for r in rows if r["visibility"] not in ("ground level", "never")]
    print(f"places: {len(rows)}  visible from the ground {len(vis)}"
          f"  rooftop-only {len(roof)}"
          f"  never {len(rows)-len(vis)-len(roof)}")
    print(f"  furthest place with a ground-level view: "
          f"{max(vis, key=lambda r: r['distance_km'])['name']} "
          f"at {max(vis, key=lambda r: r['distance_km'])['distance_km']} km")

    # ---- buildings
    bf = load_geojson(os.path.join(OSMD, "SriLanka_buildings.geojson"))
    brows = []
    for f in bf:
        lon, lat = f["geometry"]["coordinates"]
        p = f["properties"]
        hb = p.get("height_m")
        if hb is None:
            continue
        v, tier, x, y = sample_best(lon, lat)
        if v is None or v["required_height"] == -9999.0:
            continue
        need = v["required_height"]
        if need > 300 or need <= EYE:
            gains = False
        else:
            gains = hb >= need
        if not gains:
            continue
        brows.append({
            "name": p.get("name", ""), "osm_id": p.get("osm_id"),
            "lon": round(lon, 5), "lat": round(lat, 5),
            "building_height_m": hb, "height_source": p.get("height_src", ""),
            "required_height_m": round(need, 1),
            "spare_m": round(hb - need, 1),
            "distance_km": round(v["distance_km"], 1),
            "grid": tier,
        })
    brows.sort(key=lambda r: -r["spare_m"])
    fb = os.path.join(OUTD, "SriPada_buildings_visibility.csv")
    with open(fb, "w", newline="", encoding="utf-8-sig") as f:
        if brows:
            w = csv.DictWriter(f, fieldnames=list(brows[0].keys()))
            w.writeheader(); w.writerows(brows)
    print(f"buildings tested: {len(bf):,}  roof clears where ground does not: "
          f"{len(brows):,}")
    for r in brows[:12]:
        print(f"    {r['name'][:34]:34s} {r['building_height_m']:6.0f} m "
              f"needs {r['required_height_m']:6.1f} m  @{r['distance_km']:6.1f} km")

    # geojson twin for the places table
    gj = {"type": "FeatureCollection", "features": [
        {"type": "Feature",
         "geometry": {"type": "Point", "coordinates": [r["lon"], r["lat"]]},
         "properties": r} for r in rows]}
    with open(os.path.join(OUTD, "SriPada_places_visibility.geojson"), "w",
              encoding="utf-8") as f:
        json.dump(gj, f, ensure_ascii=False)
    print(f"\nwrote {fp}\n      {fb}")
