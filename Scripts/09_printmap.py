"""
Sri Pada visibility study - Step 9: print-layout map sheet.

A2 landscape composition: the island-wide "observer height needed" surface
over hillshade, with a near-field inset around the summit, legend, scale bar,
north arrow, and a title block carrying the method's headline caveats. Renders
to both PDF (vector legend/text, crisp at any zoom) and a high-res PNG.
"""
import datetime
import os
import numpy as np
from osgeo import gdal
from qgis.core import (
    QgsApplication, QgsProject, QgsRasterLayer, QgsVectorLayer,
    QgsCoordinateReferenceSystem, QgsPalettedRasterRenderer,
    QgsPrintLayout, QgsLayoutItemMap, QgsLayoutItemLegend,
    QgsLayoutItemScaleBar, QgsLayoutItemLabel, QgsLayoutItemShape,
    QgsLayoutSize, QgsUnitTypes, QgsLayoutExporter, QgsRectangle,
    QgsMarkerSymbol, QgsSingleSymbolRenderer, QgsLegendStyle,
)
from qgis.PyQt.QtGui import QColor, QFont
from qgis.PyQt.QtCore import QRectF

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEMDIR, VS, RES = (os.path.join(ROOT, d) for d in ("DEM", "Viewshed", "Results"))

QgsApplication.setPrefixPath(r"C:\Program Files\QGIS 3.40.8\apps\qgis-ltr", True)
qgs = QgsApplication([], False)
qgs.initQgis()

prj = QgsProject.instance()
prj.setCrs(QgsCoordinateReferenceSystem("EPSG:5235"))

CLASS_COLORS = [
    (0, "#ffffb2", "Visible standing on the ground"),
    (1, "#fed976", "Needs 10 m  (3 storeys)"),
    (2, "#feb24c", "Needs 30 m  (10 storeys)"),
    (3, "#fd8d3c", "Needs 60 m  (20 storeys)"),
    (4, "#f03b20", "Needs 100 m  (33 storeys)"),
    (5, "#bd0026", "Needs 150 m  (50 storeys)"),
]


def paletted(layer, entries, opacity=1.0):
    classes = [QgsPalettedRasterRenderer.Class(v, QColor(c), lbl)
               for v, c, lbl in entries]
    layer.setRenderer(QgsPalettedRasterRenderer(layer.dataProvider(), 1, classes))
    layer.renderer().setOpacity(opacity)


def add(layer):
    prj.addMapLayer(layer)
    return layer


def masked_hillshade(hillshade_path, landmask_path, out_path):
    """Punch the ocean out of a hillshade raster.

    Flat sea (elevation exactly 0 everywhere) hillshades to one uniform grey,
    which fills the raster's whole rectangular extent - burying the actual
    coastline under a solid block instead of showing the island's true shape.
    Land is 1 in the mask; land==0 is set to nodata so ocean renders
    transparent and only the real landmass paints.
    """
    if os.path.exists(out_path):
        return out_path
    hs = gdal.Open(hillshade_path)
    lm = gdal.Open(landmask_path)
    h = hs.GetRasterBand(1).ReadAsArray()
    m = lm.GetRasterBand(1).ReadAsArray()
    if m.shape != h.shape:
        lm2 = gdal.Warp("", lm, format="MEM", width=h.shape[1], height=h.shape[0],
                        outputBounds=(hs.GetGeoTransform()[0],
                                     hs.GetGeoTransform()[3] + h.shape[0]*hs.GetGeoTransform()[5],
                                     hs.GetGeoTransform()[0] + h.shape[1]*hs.GetGeoTransform()[1],
                                     hs.GetGeoTransform()[3]),
                        resampleAlg="near")
        m = lm2.GetRasterBand(1).ReadAsArray()
    out = np.where(m == 1, h, 0).astype(np.uint8)
    drv = gdal.GetDriverByName("GTiff")
    o = drv.Create(out_path, hs.RasterXSize, hs.RasterYSize, 1, gdal.GDT_Byte,
                   options=["COMPRESS=DEFLATE", "TILED=YES"])
    o.SetGeoTransform(hs.GetGeoTransform())
    o.SetProjection(hs.GetProjection())
    b = o.GetRasterBand(1)
    b.SetNoDataValue(0)
    b.WriteArray(out)
    o.FlushCache()
    o = None
    return out_path


# ---------------------------------------------------------------- layers
hill90_masked = masked_hillshade(
    os.path.join(DEMDIR, "SriLanka_hillshade_90m.tif"),
    os.path.join(VS, "SriLanka_90m_landmask.tif"),
    os.path.join(DEMDIR, "SriLanka_hillshade_90m_masked.tif"))
hill90 = add(QgsRasterLayer(hill90_masked, "Hillshade"))

cls90 = add(QgsRasterLayer(os.path.join(VS, "SriLanka_90m_class.tif"),
                           "Height needed"))
paletted(cls90, CLASS_COLORS, opacity=0.72)

hill30_masked = masked_hillshade(
    os.path.join(DEMDIR, "SriPada_hillshade_30m.tif"),
    os.path.join(VS, "SriPada_30m_landmask.tif"),
    os.path.join(DEMDIR, "SriPada_hillshade_30m_masked.tif"))
hill30 = add(QgsRasterLayer(hill30_masked, "Hillshade near"))

cls30 = add(QgsRasterLayer(os.path.join(VS, "SriPada_30m_class.tif"),
                           "Height needed near"))
paletted(cls30, CLASS_COLORS, opacity=0.72)

summit_gj = os.path.join(RES, "SriPada_summit.geojson")
summit = add(QgsVectorLayer(summit_gj, "Summit", "ogr"))
sym = QgsMarkerSymbol.createSimple({"name": "triangle", "color": "#000000",
                                    "size": "4.5", "outline_color": "#ffffff",
                                    "outline_width": "0.7"})
summit.setRenderer(QgsSingleSymbolRenderer(sym))

# ---------------------------------------------------------------- layout
layout = QgsPrintLayout(prj)
layout.initializeDefaults()
layout.setName("Sri Pada Visibility Map")
page = layout.pageCollection().pages()[0]
page.setPageSize(QgsLayoutSize(594, 420, QgsUnitTypes.LayoutMillimeters))  # A2 landscape

# ---- background panel behind the right-hand info column. Added first so
# every item that follows stacks on top of it (layout items paint in
# add-order; adding this last, as before, silently hides the inset map and
# legend beneath an opaque rectangle).
panel = QgsLayoutItemShape(layout)
panel.setShapeType(QgsLayoutItemShape.Rectangle)
panel.attemptSetSceneRect(QRectF(440, 4, 150, 412))
panel.symbol().setColor(QColor("#18232e"))
layout.addLayoutItem(panel)

# ---- main map (island-wide)
mmain = QgsLayoutItemMap(layout)
mmain.attemptSetSceneRect(QRectF(8, 8, 430, 404))
mmain.setLayers([cls90, hill90])
# QgsLayoutItemMap.setExtent adjusts the requested extent to match the frame's
# aspect ratio, which - if the requested box has a different aspect - crops
# rather than pads. Pre-shape the box to the frame's own aspect (430:404) so
# the island isn't clipped at the north or south ends.
frame_aspect = 430.0 / 404.0
isl_cx, isl_cy = 497_500.0, 592_500.0     # centre of the true DEM footprint
isl_h = 620_000.0                         # generous N-S margin
isl_w = isl_h * frame_aspect
req = QgsRectangle(isl_cx - isl_w / 2, isl_cy - isl_h / 2,
                   isl_cx + isl_w / 2, isl_cy + isl_h / 2)
print(f"requested extent: {req.toString()}  aspect={req.width()/req.height():.4f}")
mmain.setExtent(req)
mmain.setFrameEnabled(True)
layout.addLayoutItem(mmain)
got = mmain.extent()
print(f"actual extent:    {got.toString()}  aspect={got.width()/got.height():.4f}")
print(f"class raster extent: {cls90.extent().toString()}")

# ---- inset map (near field)
minset = QgsLayoutItemMap(layout)
minset.attemptSetSceneRect(QRectF(444, 8, 142, 142))
minset.setLayers([cls30, hill30, summit])
near_ext = QgsRectangle(469682.5 - 40000, 478848.2 - 40000,
                        469682.5 + 40000, 478848.2 + 40000)
minset.setExtent(near_ext)
minset.setFrameEnabled(True)
minset.setBackgroundColor(QColor("#0f1720"))
layout.addLayoutItem(minset)

lbl_inset = QgsLayoutItemLabel(layout)
lbl_inset.setText("NEAR FIELD \u2014 80 km around the summit")
lbl_inset.setFont(QFont("Arial", 8, QFont.Bold))
lbl_inset.setFontColor(QColor("#e8eef4"))
lbl_inset.attemptSetSceneRect(QRectF(444, 152, 142, 6))
layout.addLayoutItem(lbl_inset)

# ---- legend
leg = QgsLayoutItemLegend(layout)
leg.setLinkedMap(mmain)
leg.setTitle("Observer height needed to see Sri Pada")
leg.setAutoUpdateModel(False)
root = leg.model().rootGroup()
root.clear()
leg_node = root.addLayer(cls90)
leg.setStyleFont(QgsLegendStyle.Title, QFont("Arial", 10, QFont.Bold))
leg.setStyleFont(QgsLegendStyle.Subgroup, QFont("Arial", 9))
leg.setStyleFont(QgsLegendStyle.SymbolLabel, QFont("Arial", 9))
# default legend text is near-black, unreadable against the dark info panel
leg.setFontColor(QColor("#e8eef4"))
# QGIS raster legends insert a "Band 1 (Gray)" sub-header above the class
# swatches by default; blank it out rather than confuse the reader with an
# irrelevant internal band label.
for n in leg.model().layerOriginalLegendNodes(leg_node):
    if "Band 1" in n.data(0):
        n.setUserLabel(" ")
leg.attemptSetSceneRect(QRectF(444, 162, 142, 92))
leg.setBackgroundColor(QColor(24, 35, 46, 235))
layout.addLayoutItem(leg)

# ---- scale bar
sb = QgsLayoutItemScaleBar(layout)
sb.setLinkedMap(mmain)
sb.setStyle("Single Box")
sb.setUnits(QgsUnitTypes.DistanceKilometers)
sb.setNumberOfSegments(4)
sb.setNumberOfSegmentsLeft(0)
sb.setUnitsPerSegment(50)
sb.setUnitLabel("km")
sb.setFont(QFont("Arial", 8))
sb.update()
sb.attemptSetSceneRect(QRectF(20, 392, 90, 14))
layout.addLayoutItem(sb)

# ---- north arrow (label glyph, avoids external svg dependency). Sits over
# the map frame, which is light hillshade grey / white margin, not the dark
# info panel, so it needs a dark colour to actually be visible.
na = QgsLayoutItemLabel(layout)
na.setText("N\n\u25B2")
na.setFont(QFont("Arial", 11, QFont.Bold))
na.setFontColor(QColor("#1a1a1a"))
na.attemptSetSceneRect(QRectF(400, 14, 16, 20))
layout.addLayoutItem(na)


# ---- title block
def label(text, x, y, w, h, size=10, bold=False, color="#1a1a1a"):
    it = QgsLayoutItemLabel(layout)
    it.setText(text)
    it.setFont(QFont("Arial", size, QFont.Bold if bold else QFont.Normal))
    it.setFontColor(QColor(color))
    it.attemptSetSceneRect(QRectF(x, y, w, h))
    layout.addLayoutItem(it)
    return it


label("SRI PADA VISIBILITY", 444, 262, 142, 10, size=15, bold=True,
      color="#e8eef4")
label("Adam's Peak  \u00b7  6.8094N 80.4994E  \u00b7  2192 m DSM / 2243 m surveyed",
      444, 273, 142, 8, size=8, color="#9fb0c0")

caveat = (
    "Line of sight from a 30 m / 90 m Copernicus GLO-30 DSM, corrected for "
    "earth curvature and atmospheric refraction (k=0.13). Canopy is included "
    "in the surface: it blocks realistically but also lifts forested "
    "observers. Colour shows the minimum observer height above ground needed "
    "to clear the terrain between that point and the summit; the same "
    "statistic run in reverse is a tall-building viewshed. Haze is not "
    "modelled: treat anything beyond about 100 km as geometric only, not a "
    "real sighting. CRS: EPSG:5235 (SLD99)."
)
label(caveat, 444, 284, 142, 60, size=7, color="#c7d3dd")

label("15.4% of Sri Lanka's land sees the summit from the ground; 37.2% "
      "with 150 m of height beneath it; 62.8% never, at any height.",
      444, 346, 142, 20, size=8, bold=True, color="#e8eef4")

label("Compiled " + datetime.date.today().isoformat(), 444, 406, 142, 6,
      size=7, color="#7a8b9c")

prj.layoutManager().addLayout(layout)

# ---------------------------------------------------------------- export
out_pdf = os.path.join(ROOT, "SriPada_Visibility_Map_A2.pdf")
out_png = os.path.join(ROOT, "SriPada_Visibility_Map_A2.png")

exp = QgsLayoutExporter(layout)
r1 = exp.exportToPdf(out_pdf, QgsLayoutExporter.PdfExportSettings())
print("PDF export result:", r1, "->", out_pdf)

img_settings = QgsLayoutExporter.ImageExportSettings()
img_settings.dpi = 200
r2 = exp.exportToImage(out_png, img_settings)
print("PNG export result:", r2, "->", out_png)

qgs.exitQgis()
