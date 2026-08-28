"""
Sri Pada visibility study - Step 10: single 30 m island-wide base.

Replaces the old two-tier (30 m near field / 90 m island) arrangement. COP30's
native posting is 30 m, so this is the finest the source supports; nothing in
the pipeline resamples coarser than the data from here on.

Outputs, all EPSG:5235 at 30 m:
  SriLanka_DEM_30m.tif        9501 x 15834, the working surface
  SriLanka_landmask_30m.tif   land = 1 (Copernicus codes sea as exact 0.0)
  SriLanka_hillshade_30m.tif  full hillshade
  SriLanka_hillshade_30m_masked.tif   ocean punched out, for map rendering
"""
import os
import numpy as np
from osgeo import gdal

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEM = os.path.join(ROOT, "DEM")
VRT = os.path.join(DEM, "SriLanka_COP30_wgs84.vrt")

ISLAND = (355_000, 354_980, 640_030, 830_000)   # xmin ymin xmax ymax
CO = ["COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=YES",
      "NUM_THREADS=ALL_CPUS"]

dem30 = os.path.join(DEM, "SriLanka_DEM_30m.tif")
if not os.path.exists(dem30):
    print("warping island to 30 m ...", flush=True)
    gdal.Warp(dem30, VRT, dstSRS="EPSG:5235", xRes=30, yRes=30,
              outputBounds=ISLAND, resampleAlg="bilinear",
              outputType=gdal.GDT_Float32, creationOptions=CO,
              multithread=True, warpMemoryLimit=512)
ds = gdal.Open(dem30)
print(f"DEM {ds.RasterXSize} x {ds.RasterYSize} @30 m "
      f"({ds.RasterXSize*ds.RasterYSize/1e6:.0f} M cells)")

# ---- land mask: nearest-neighbour so the exact sea zeros survive
lm30 = os.path.join(DEM, "SriLanka_landmask_30m.tif")
if not os.path.exists(lm30):
    print("land mask ...", flush=True)
    tmp = gdal.Warp("", VRT, format="VRT", dstSRS="EPSG:5235",
                    xRes=30, yRes=30, outputBounds=ISLAND, resampleAlg="near",
                    outputType=gdal.GDT_Float32)
    drv = gdal.GetDriverByName("GTiff")
    o = drv.Create(lm30, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte,
                   options=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES"])
    o.SetGeoTransform(ds.GetGeoTransform())
    o.SetProjection(ds.GetProjection())
    band = o.GetRasterBand(1)
    nland = 0
    BL = 2048
    for r0 in range(0, ds.RasterYSize, BL):
        h = min(BL, ds.RasterYSize - r0)
        z = tmp.GetRasterBand(1).ReadAsArray(0, r0, ds.RasterXSize, h)
        m = (z != 0.0).astype(np.uint8)
        band.WriteArray(m, 0, r0)
        nland += int(m.sum())
    band.FlushCache()
    o = None
    print(f"  land {nland*900/1e6:,.0f} km2   (surveyed 65,610 km2)")

# ---- hillshades
hs = os.path.join(DEM, "SriLanka_hillshade_30m.tif")
if not os.path.exists(hs):
    print("hillshade ...", flush=True)
    gdal.DEMProcessing(hs, dem30, "hillshade", azimuth=315, altitude=45,
                       zFactor=1.5, computeEdges=True,
                       creationOptions=["COMPRESS=DEFLATE", "TILED=YES",
                                        "BIGTIFF=YES"])

hsm = os.path.join(DEM, "SriLanka_hillshade_30m_masked.tif")
if not os.path.exists(hsm):
    print("masked hillshade ...", flush=True)
    h_ds = gdal.Open(hs)
    l_ds = gdal.Open(lm30)
    drv = gdal.GetDriverByName("GTiff")
    o = drv.Create(hsm, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte,
                   options=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES"])
    o.SetGeoTransform(ds.GetGeoTransform())
    o.SetProjection(ds.GetProjection())
    b = o.GetRasterBand(1)
    b.SetNoDataValue(0)
    BL = 2048
    for r0 in range(0, ds.RasterYSize, BL):
        hgt = min(BL, ds.RasterYSize - r0)
        hh = h_ds.GetRasterBand(1).ReadAsArray(0, r0, ds.RasterXSize, hgt)
        ll = l_ds.GetRasterBand(1).ReadAsArray(0, r0, ds.RasterXSize, hgt)
        b.WriteArray(np.where(ll == 1, hh, 0).astype(np.uint8), 0, r0)
    b.FlushCache()
    o = None

print("done")
