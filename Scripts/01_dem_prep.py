"""
Sri Pada visibility study - Step 1: DEM preparation.

Mosaics the 12 Copernicus GLO-30 tiles covering Sri Lanka and produces two
working grids in EPSG:5235 (SLD99 / Sri Lanka Grid 1999, metres):

  SriLanka_DEM_90m.tif   island-wide, 90 m  -> far-field viewshed
  SriPada_DEM_30m.tif    110 km around the summit, 30 m -> near-field detail

Ocean / missing-tile areas resolve to 0 m (sea level), which is the correct
surface for a line of sight crossing water.
"""
import os
from osgeo import gdal, osr

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEM = os.path.join(ROOT, "DEM")

SUMMIT_LON, SUMMIT_LAT = 80.499444, 6.809444   # Sri Pada / Adam's Peak
NEARFIELD_RADIUS_M = 110_000                    # 30 m tier half-width

# ---------------------------------------------------------------- summit XY
wgs = osr.SpatialReference(); wgs.ImportFromEPSG(4326)
wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
sld = osr.SpatialReference(); sld.ImportFromEPSG(5235)
sld.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
to_sld = osr.CoordinateTransformation(wgs, sld)
sx, sy, _ = to_sld.TransformPoint(SUMMIT_LON, SUMMIT_LAT)
print(f"Summit EPSG:5235  X={sx:.1f}  Y={sy:.1f}")

# ------------------------------------------------------------------ mosaic
tiles = sorted(os.path.join(DEM, f) for f in os.listdir(DEM)
               if f.startswith("COP30_") and f.endswith(".tif"))
print(f"{len(tiles)} tiles")
vrt = os.path.join(DEM, "SriLanka_COP30_wgs84.vrt")
gdal.BuildVRT(vrt, tiles)

# Island bounds in EPSG:5235, generous margin so nothing is clipped.
ISLAND = (355_000, 355_000, 640_000, 830_000)   # xmin ymin xmax ymax

co = ["COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"]

out90 = os.path.join(DEM, "SriLanka_DEM_90m.tif")
gdal.Warp(out90, vrt, dstSRS="EPSG:5235", xRes=90, yRes=90,
          outputBounds=ISLAND, resampleAlg="average",
          outputType=gdal.GDT_Float32, dstNodata=None,
          creationOptions=co, multithread=True)
print("wrote", out90)

nf = (sx - NEARFIELD_RADIUS_M, sy - NEARFIELD_RADIUS_M,
      sx + NEARFIELD_RADIUS_M, sy + NEARFIELD_RADIUS_M)
out30 = os.path.join(DEM, "SriPada_DEM_30m.tif")
gdal.Warp(out30, vrt, dstSRS="EPSG:5235", xRes=30, yRes=30,
          outputBounds=nf, resampleAlg="bilinear",
          outputType=gdal.GDT_Float32, dstNodata=None,
          creationOptions=co, multithread=True)
print("wrote", out30)

# --------------------------------------------------------------- summarise
for p in (out90, out30):
    ds = gdal.Open(p)
    b = ds.GetRasterBand(1)
    mn, mx, mean, sd = b.ComputeStatistics(False)
    gt = ds.GetGeoTransform()
    print(f"{os.path.basename(p):22s} {ds.RasterXSize}x{ds.RasterYSize} "
          f"@{gt[1]:g} m   min={mn:.1f} max={mx:.1f} mean={mean:.1f}")
    ds = None

# summit elevation as the DEM sees it
ds = gdal.Open(out30)
gt = ds.GetGeoTransform()
col = int((sx - gt[0]) / gt[1]); row = int((sy - gt[3]) / gt[5])
win = ds.GetRasterBand(1).ReadAsArray(col - 5, row - 5, 11, 11)
print(f"DEM elevation at summit cell = {win[5,5]:.1f} m ; "
      f"max in 11x11 window = {win.max():.1f} m")
