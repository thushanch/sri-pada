"""
Sri Pada visibility study - Step 12: interactive map + "view from here".

Builds a single self-contained HTML with two tabs:

  MAP    the 30 m visibility surface over a choice of basemaps, including an
         offline one baked from the study's own hillshade so the page still
         works with no network at all.

  VIEW   click anywhere and the page renders the skyline you would actually
         see looking toward Sri Pada from that spot: intervening ridges drawn
         far-to-near with haze, the summit marked, and a live eye-height
         slider so the tall-building result can be felt rather than read.

Terrain for the panorama is embedded as PNG-packed height grids (elevation in
the R/G channels, 0.1 m steps). PNG's DEFLATE squeezes the flat ocean to
almost nothing, which is why this is far smaller than raw base64 arrays.

The visibility VERDICT always comes from the full 30 m analysis raster; the
panorama's ridge SHAPE is drawn from a ~220 m resample, which is stated in the
UI so the two are never confused.
"""
import base64
import json
import os
import numpy as np
from osgeo import gdal

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEMDIR, VS, SOL, RES = (os.path.join(ROOT, d)
                        for d in ("DEM", "Viewshed", "Solar", "Results"))
WEB = os.path.join(ROOT, "Web")
os.makedirs(WEB, exist_ok=True)

SUMMIT_LON, SUMMIT_LAT, SUMMIT_Z = 80.499444, 6.809444, 2192.0
GRID_DEG = 0.002              # ~222 m, panorama terrain
MAXDIM = 2600

CLASS_RGBA = {0: (255, 255, 178, 205), 1: (254, 217, 118, 205),
              2: (254, 178, 76, 205), 3: (253, 141, 60, 205),
              4: (240, 59, 32, 205), 5: (189, 0, 38, 205)}


def png_bytes(arr_rgba):
    h, w, nb = arr_rgba.shape
    mem = gdal.GetDriverByName("MEM").Create("", w, h, nb, gdal.GDT_Byte)
    for i in range(nb):
        mem.GetRasterBand(i + 1).WriteArray(arr_rgba[:, :, i])
    tmp = f"/vsimem/p_{id(arr_rgba)}.png"
    gdal.GetDriverByName("PNG").CreateCopy(tmp, mem, options=["ZLEVEL=9"])
    f = gdal.VSIFOpenL(tmp, "rb")
    gdal.VSIFSeekL(f, 0, 2)
    n = gdal.VSIFTellL(f)
    gdal.VSIFSeekL(f, 0, 0)
    data = gdal.VSIFReadL(1, n, f)
    gdal.VSIFCloseL(f)
    gdal.Unlink(tmp)
    return data


def packed_grid(path, scale=10.0, res=GRID_DEG, alg="bilinear"):
    """Warp to lat/lon and pack values into PNG R/G as uint16."""
    w = gdal.Warp("", gdal.Open(path), format="MEM", dstSRS="EPSG:4326",
                  xRes=res, yRes=res, resampleAlg=alg)
    a = w.GetRasterBand(1).ReadAsArray().astype(np.float64)
    nd = w.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        a = np.where(a == nd, -3276.0, a)
    a = np.where(np.isfinite(a), a, -3276.0)
    v = np.clip(np.round(a * scale), -32768, 32767).astype(np.int32) + 32768
    rgb = np.zeros((a.shape[0], a.shape[1], 3), np.uint8)
    rgb[:, :, 0] = (v >> 8) & 0xFF
    rgb[:, :, 1] = v & 0xFF
    gt = w.GetGeoTransform()
    data = png_bytes(rgb)
    print(f"  {os.path.basename(path):42s} {a.shape[1]}x{a.shape[0]}  "
          f"{len(data)/1e6:5.2f} MB")
    return {"b64": base64.b64encode(data).decode(),
            "x0": gt[0], "y0": gt[3], "dx": gt[1], "dy": gt[5],
            "scale": scale}


def overlay(path, lut, alg="mode"):
    ds = gdal.Open(path)
    w = gdal.Warp("", ds, format="MEM", dstSRS="EPSG:3857", resampleAlg=alg)
    sc = max(w.RasterXSize, w.RasterYSize) / MAXDIM
    if sc > 1:
        w = gdal.Warp("", w, format="MEM", width=int(w.RasterXSize / sc),
                      height=int(w.RasterYSize / sc), resampleAlg=alg)
    a = w.GetRasterBand(1).ReadAsArray()
    rgba = np.zeros((a.shape[0], a.shape[1], 4), np.uint8)
    for v, c in lut.items():
        rgba[a == v] = c
    gt = w.GetGeoTransform()
    from osgeo import osr
    m = osr.SpatialReference(); m.ImportFromEPSG(3857)
    m.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    g = osr.SpatialReference(); g.ImportFromEPSG(4326)
    g.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(m, g)
    x1 = gt[0] + gt[1] * w.RasterXSize
    y1 = gt[3] + gt[5] * w.RasterYSize
    (lo1, la1, _) = tr.TransformPoint(gt[0], y1)
    (lo2, la2, _) = tr.TransformPoint(x1, gt[3])
    data = png_bytes(rgba)
    print(f"  {os.path.basename(path):42s} {a.shape[1]}x{a.shape[0]}  "
          f"{len(data)/1e6:5.2f} MB")
    return {"b64": base64.b64encode(data).decode(),
            "bounds": [[la1, lo1], [la2, lo2]]}


def hillshade_overlay(path, landmask):
    """Greyscale relief as an RGBA overlay, ocean transparent - the offline
    basemap, so the page is useful with no network."""
    ds = gdal.Open(path)
    w = gdal.Warp("", ds, format="MEM", dstSRS="EPSG:3857",
                  resampleAlg="bilinear")
    sc = max(w.RasterXSize, w.RasterYSize) / MAXDIM
    if sc > 1:
        w = gdal.Warp("", w, format="MEM", width=int(w.RasterXSize / sc),
                      height=int(w.RasterYSize / sc), resampleAlg="bilinear")
    a = w.GetRasterBand(1).ReadAsArray()
    wgt = w.GetGeoTransform()
    # keep a reference: chaining off gdal.Warp() frees the dataset before the
    # band is read, which surfaces as a confusing SWIG type error
    lm_ds = gdal.Warp("", gdal.Open(landmask), format="MEM", dstSRS="EPSG:3857",
                      width=w.RasterXSize, height=w.RasterYSize,
                      outputBounds=[wgt[0],
                                    wgt[3] + wgt[5] * w.RasterYSize,
                                    wgt[0] + wgt[1] * w.RasterXSize,
                                    wgt[3]],
                      resampleAlg="near")
    lm = lm_ds.GetRasterBand(1).ReadAsArray()
    rgba = np.zeros((a.shape[0], a.shape[1], 4), np.uint8)
    tone = (140 + (a.astype(np.int32) - 140) * 0.55).clip(60, 245).astype(np.uint8)
    rgba[:, :, 0] = tone
    rgba[:, :, 1] = tone
    rgba[:, :, 2] = (tone.astype(np.int32) * 0.97).clip(0, 255).astype(np.uint8)
    rgba[:, :, 3] = np.where(lm == 1, 255, 0).astype(np.uint8)
    gt = w.GetGeoTransform()
    from osgeo import osr
    m = osr.SpatialReference(); m.ImportFromEPSG(3857)
    m.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    g = osr.SpatialReference(); g.ImportFromEPSG(4326)
    g.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    tr = osr.CoordinateTransformation(m, g)
    x1 = gt[0] + gt[1] * w.RasterXSize
    y1 = gt[3] + gt[5] * w.RasterYSize
    (lo1, la1, _) = tr.TransformPoint(gt[0], y1)
    (lo2, la2, _) = tr.TransformPoint(x1, gt[3])
    data = png_bytes(rgba)
    print(f"  offline relief basemap                     "
          f"{a.shape[1]}x{a.shape[0]}  {len(data)/1e6:5.2f} MB")
    return {"b64": base64.b64encode(data).decode(),
            "bounds": [[la1, lo1], [la2, lo2]]}


def declination_table():
    src = open(os.path.join(ROOT, "Scripts", "04_solar.py"),
               encoding="utf-8").read().split('if __name__ == "__main__":')[0]
    ns = {}
    exec(compile(src, "s4", "exec"), ns)
    _, dec = ns["declination_series"](2026)
    return [round(float(d), 4) for d in dec]


if __name__ == "__main__":
    print("overlays ...")
    layers = {}
    layers["height_needed"] = overlay(
        os.path.join(VS, "SriLanka_30m_class.tif"), CLASS_RGBA)
    relief = hillshade_overlay(
        os.path.join(DEMDIR, "SriLanka_hillshade_30m_masked.tif"),
        os.path.join(DEMDIR, "SriLanka_landmask_30m.tif"))

    print("terrain grids ...")
    grids = {
        "z": packed_grid(os.path.join(DEMDIR, "SriLanka_DEM_30m.tif"),
                         scale=10.0),
        "h": packed_grid(os.path.join(VS, "SriLanka_30m_required_height.tif"),
                         scale=1.0),
    }

    # Every named settlement, with the answer already attached, so the town
    # search can colour each hit before the user even clicks it.
    places = []
    pv = os.path.join(RES, "SriPada_places_visibility.geojson")
    if os.path.exists(pv):
        with open(pv, encoding="utf-8") as f:
            for ft in json.load(f)["features"]:
                p = ft["properties"]
                if not p.get("name"):
                    continue
                places.append({"n": p["name"], "lo": p["lon"], "la": p["lat"],
                               "d": p["distance_km"], "v": p["visibility"]})
    rank = {"city": 0, "town": 1, "suburb": 2, "village": 3}
    places.sort(key=lambda r: r["d"])
    print(f"  {len(places)} named places embedded for search")

    brand = None
    bp = os.path.join(ROOT, "brand.json")
    if os.path.exists(bp):
        brand = json.load(open(bp, encoding="utf-8"))
        print("  brand assets embedded (wordmark + mark)")
    else:
        print("  ! brand.json missing - run SriLankaPeaks/Scripts/24_branding.py")

    payload = {"layers": layers, "relief": relief, "grids": grids,
               "brand": brand,
               "places": places, "dec": declination_table(),
               "summit": {"lon": SUMMIT_LON, "lat": SUMMIT_LAT, "z": SUMMIT_Z}}

    sd = os.path.join(ROOT, "Scripts")
    tpl = open(os.path.join(sd, "web_template30.html"), encoding="utf-8").read()

    # Inline Leaflet so the page needs no network at all: with the library on
    # a CDN, "offline basemap" only meant the tiles - the map itself would
    # still fail to load.
    for token, fn in (("/*__LEAFLET_CSS__*/", "leaflet.css"),
                      ("/*__LEAFLET_JS__*/", "leaflet.js")):
        p = os.path.join(sd, fn)
        if not os.path.exists(p):
            raise SystemExit(f"missing {p} - fetch it from unpkg first")
        tpl = tpl.replace(token, open(p, encoding="utf-8").read())

    out = os.path.join(WEB, "SriPada_Visibility_Map.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(tpl.replace("/*__DATA__*/", json.dumps(payload)))
    print(f"\nwrote {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
