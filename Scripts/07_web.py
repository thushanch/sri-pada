"""
Sri Pada visibility study - Step 7: self-contained interactive web map.

Raster layers are warped to Web Mercator, colourised, and embedded as base64
PNGs, so the file works from disk with no server and no sidecar assets (it does
fetch OSM basemap tiles when online).

Two coarse Int16 grids are embedded as well - required observer height and
ground elevation - which is enough for the page to answer, for any point the
user clicks: distance, bearing, whether the peak is visible, how high you would
have to be, how large it looks, and the two dates the rising sun comes up from
behind it. The sunrise dates are solved live in the browser from the same
declination identity the raster pipeline uses, so the map and the click
readout cannot drift apart.
"""
import base64
import json
import os
import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEMDIR, VS, SOL, RES = (os.path.join(ROOT, d)
                        for d in ("DEM", "Viewshed", "Solar", "Results"))
WEB = os.path.join(ROOT, "Web")
os.makedirs(WEB, exist_ok=True)

SUMMIT_LON, SUMMIT_LAT, SUMMIT_Z = 80.499444, 6.809444, 2192.0
MAXDIM = 2400

CLASS_RGBA = {0: (255, 255, 178, 205), 1: (254, 217, 118, 205),
              2: (254, 178, 76, 205), 3: (253, 141, 60, 205),
              4: (240, 59, 32, 205), 5: (189, 0, 38, 205)}
GRADE_RGBA = {1: (8, 48, 107, 215), 2: (33, 113, 181, 205),
              3: (66, 146, 198, 195), 4: (158, 202, 225, 185),
              5: (222, 235, 247, 175)}
SHADOW_RGBA = {1: (84, 39, 143, 120), 3: (253, 174, 107, 230)}


def warp3857(path, max_dim=MAXDIM, alg="near"):
    ds = gdal.Open(path)
    w = gdal.Warp("", ds, format="MEM", dstSRS="EPSG:3857", resampleAlg=alg)
    sc = max(w.RasterXSize, w.RasterYSize) / max_dim
    if sc > 1:
        w = gdal.Warp("", w, format="MEM",
                      width=int(w.RasterXSize / sc),
                      height=int(w.RasterYSize / sc), resampleAlg=alg)
    return w


def bounds_latlon(ds):
    gt = ds.GetGeoTransform()
    xs = [gt[0], gt[0] + gt[1] * ds.RasterXSize]
    ys = [gt[3] + gt[5] * ds.RasterYSize, gt[3]]
    m = osr.SpatialReference(); m.ImportFromEPSG(3857)
    m.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    g = osr.SpatialReference(); g.ImportFromEPSG(4326)
    g.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(m, g)
    (lo1, la1, _), (lo2, la2, _) = (tr.TransformPoint(xs[0], ys[0]),
                                    tr.TransformPoint(xs[1], ys[1]))
    return [[la1, lo1], [la2, lo2]]


def png_b64(rgba):
    h, w, _ = rgba.shape
    mem = gdal.GetDriverByName("MEM").Create("", w, h, 4, gdal.GDT_Byte)
    for i in range(4):
        mem.GetRasterBand(i + 1).WriteArray(rgba[:, :, i])
    tmp = f"/vsimem/o_{id(rgba)}.png"
    gdal.GetDriverByName("PNG").CreateCopy(tmp, mem, options=["ZLEVEL=9"])
    f = gdal.VSIFOpenL(tmp, "rb")
    gdal.VSIFSeekL(f, 0, 2); n = gdal.VSIFTellL(f); gdal.VSIFSeekL(f, 0, 0)
    data = gdal.VSIFReadL(1, n, f)
    gdal.VSIFCloseL(f); gdal.Unlink(tmp)
    return base64.b64encode(data).decode(), n


def categorical_layer(path, lut, alg="mode"):
    ds = warp3857(path, alg=alg)
    a = ds.GetRasterBand(1).ReadAsArray()
    rgba = np.zeros((a.shape[0], a.shape[1], 4), np.uint8)
    for v, (r, g, b, al) in lut.items():
        m = a == v
        rgba[m] = (r, g, b, al)
    b64, n = png_b64(rgba)
    return {"b64": b64, "bounds": bounds_latlon(ds)}, n


def ramp_layer(path, vmin, vmax, stops, alg="bilinear"):
    ds = warp3857(path, alg=alg)
    a = ds.GetRasterBand(1).ReadAsArray().astype(np.float32)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    valid = np.isfinite(a) & (a != (nd if nd is not None else -9999.0))
    t = np.clip((a - vmin) / (vmax - vmin), 0, 1)
    xs = np.array([s[0] for s in stops], np.float32)
    cols = np.array([s[1] for s in stops], np.float32)
    rgba = np.zeros((a.shape[0], a.shape[1], 4), np.uint8)
    for i in range(3):
        rgba[:, :, i] = np.interp(t, xs, cols[:, i]).astype(np.uint8)
    rgba[:, :, 3] = np.where(valid, 210, 0).astype(np.uint8)
    b64, n = png_b64(rgba)
    return {"b64": b64, "bounds": bounds_latlon(ds)}, n


def query_grid(path, res_deg=0.0045, as_int=True, alg="bilinear"):
    ds = gdal.Open(path)
    w = gdal.Warp("", ds, format="MEM", dstSRS="EPSG:4326",
                  xRes=res_deg, yRes=res_deg, resampleAlg=alg)
    a = w.GetRasterBand(1).ReadAsArray().astype(np.float32)
    a = np.where(np.isfinite(a), a, -9999.0)
    a = np.clip(a, -32000, 32000).astype(np.int16)
    gt = w.GetGeoTransform()
    return {
        "b64": base64.b64encode(a.tobytes()).decode(),
        "w": int(w.RasterXSize), "h": int(w.RasterYSize),
        "x0": gt[0], "y0": gt[3], "dx": gt[1], "dy": gt[5],
    }


def declination_table():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "s4", os.path.join(ROOT, "Scripts", "04_solar.py"))
    src = open(spec.origin, encoding="utf-8").read().split(
        'if __name__ == "__main__":')[0]
    ns = {}
    exec(compile(src, "s4", "exec"), ns)
    doy, dec = ns["declination_series"](2026)
    return [round(float(d), 4) for d in dec]


if __name__ == "__main__":
    layers, total = {}, 0
    print("rendering overlays ...")

    L, n = categorical_layer(os.path.join(VS, "SriLanka_90m_class.tif"),
                             CLASS_RGBA)
    layers["height_needed"] = L; total += n
    print(f"  height needed          {n/1e6:5.2f} MB")

    L, n = categorical_layer(os.path.join(VS, "SriLanka_90m_grade.tif"),
                             GRADE_RGBA)
    layers["grade"] = L; total += n
    print(f"  how strongly it reads  {n/1e6:5.2f} MB")

    p = os.path.join(SOL, "SriPada_30m_align_doy1.tif")
    if os.path.exists(p):
        L, n = ramp_layer(p, 1, 365,
                          [(0.0, (69, 117, 180)), (0.22, (145, 191, 219)),
                           (0.47, (254, 224, 144)), (0.73, (252, 141, 89)),
                           (1.0, (215, 48, 39))])
        layers["align"] = L; total += n
        print(f"  sunrise alignment date {n/1e6:5.2f} MB")

    shp = sorted(f for f in os.listdir(SOL)
                 if f.startswith("SriPada_shadow_") and f.endswith("alt0.5.tif"))
    for f in shp:
        key = "shadow_" + f.split("_")[2]
        L, n = categorical_layer(os.path.join(SOL, f), SHADOW_RGBA)
        layers[key] = L; total += n
        print(f"  {f[:-4]:34s} {n/1e6:5.2f} MB")

    print("embedding query grids ...")
    q_h = query_grid(os.path.join(VS, "SriLanka_90m_required_height.tif"))
    q_z = query_grid(os.path.join(DEMDIR, "SriLanka_DEM_90m.tif"))
    print(f"  grid {q_h['w']}x{q_h['h']}")

    places = []
    pv = os.path.join(RES, "SriPada_places_visibility.geojson")
    if os.path.exists(pv):
        with open(pv, encoding="utf-8") as f:
            for ft in json.load(f)["features"]:
                p = ft["properties"]
                if p.get("place") in ("city", "town") or \
                        p.get("visibility") == "ground level":
                    places.append({
                        "n": p.get("name", ""), "lo": p["lon"], "la": p["lat"],
                        "d": p["distance_km"], "v": p["visibility"],
                        "h": p["required_height_m"],
                        "a": p.get("angular_height_deg", ""),
                        "s1": p.get("sunrise_align_1", ""),
                        "s2": p.get("sunrise_align_2", ""),
                    })
    places = sorted(places, key=lambda r: r["d"])[:900]
    print(f"  {len(places)} places embedded")

    payload = {
        "layers": layers, "qh": q_h, "qz": q_z, "places": places,
        "dec": declination_table(),
        "summit": {"lon": SUMMIT_LON, "lat": SUMMIT_LAT, "z": SUMMIT_Z},
    }

    tpl = open(os.path.join(ROOT, "Scripts", "web_template.html"),
               encoding="utf-8").read()
    html = tpl.replace("/*__DATA__*/", json.dumps(payload))
    out = os.path.join(WEB, "SriPada_Visibility_Map.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\nwrote {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
