"""
Sri Pada visibility study - Step 6: hillshades and the QGIS project.

Builds Sri Pada Visibility.qgz with the layers stacked and styled, an OSM
basemap underneath and the summit marked. Run headless; open the .qgz in
QGIS 3.40 afterwards.
"""
import json
import os
from osgeo import gdal

gdal.UseExceptions()

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
DEMDIR, VS, OSMD, SOL = (os.path.join(ROOT, d)
                         for d in ("DEM", "Viewshed", "OSM", "Solar"))
RES = os.path.join(ROOT, "Results")

SUMMIT_LON, SUMMIT_LAT = 80.499444, 6.809444

CLASS_COLORS = [
    (0, "#ffffb2", "Visible standing on the ground"),
    (1, "#fed976", "Needs 10 m (3 storeys)"),
    (2, "#feb24c", "Needs 30 m (10 storeys)"),
    (3, "#fd8d3c", "Needs 60 m (20 storeys)"),
    (4, "#f03b20", "Needs 100 m (33 storeys)"),
    (5, "#bd0026", "Needs 150 m (50 storeys)"),
]
GRADE_COLORS = [
    (1, "#08306b", "Prominent  >= 2 deg"),
    (2, "#2171b5", "Clear  1-2 deg"),
    (3, "#4292c6", "Distinct  0.5-1 deg"),
    (4, "#9ecae1", "Faint  0.2-0.5 deg"),
    (5, "#deebf7", "Marginal  < 0.2 deg"),
]


def hillshades():
    for src, dst in ((("SriPada_DEM_30m.tif"), "SriPada_hillshade_30m.tif"),
                     (("SriLanka_DEM_90m.tif"), "SriLanka_hillshade_90m.tif")):
        p = os.path.join(DEMDIR, dst)
        if os.path.exists(p):
            print("  hillshade exists:", dst)
            continue
        gdal.DEMProcessing(p, os.path.join(DEMDIR, src), "hillshade",
                           azimuth=315, altitude=45, zFactor=1.5,
                           computeEdges=True,
                           creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])
        print("  wrote", dst)


def summit_geojson():
    p = os.path.join(RES, "SriPada_summit.geojson")
    gj = {"type": "FeatureCollection", "features": [{
        "type": "Feature",
        "geometry": {"type": "Point",
                     "coordinates": [SUMMIT_LON, SUMMIT_LAT]},
        "properties": {"name": "Sri Pada / Adam's Peak",
                       "dem_elev_m": 2192.0, "surveyed_elev_m": 2243.0}}]}
    with open(p, "w", encoding="utf-8") as f:
        json.dump(gj, f)
    return p


def build_project():
    from qgis.core import (QgsApplication, QgsProject, QgsRasterLayer,
                           QgsVectorLayer, QgsCoordinateReferenceSystem,
                           QgsPalettedRasterRenderer, QgsRasterShader,
                           QgsColorRampShader, QgsSingleBandPseudoColorRenderer,
                           QgsMarkerSymbol, QgsSingleSymbolRenderer,
                           QgsLayerTreeGroup)
    from qgis.PyQt.QtGui import QColor

    QgsApplication.setPrefixPath(r"C:\Program Files\QGIS 3.40.8\apps\qgis-ltr",
                                 True)
    qgs = QgsApplication([], False)
    qgs.initQgis()

    prj = QgsProject.instance()
    prj.setCrs(QgsCoordinateReferenceSystem("EPSG:5235"))
    root = prj.layerTreeRoot()

    def add(layer, group=None, visible=True):
        prj.addMapLayer(layer, False)
        node = (group or root).addLayer(layer)
        node.setItemVisibilityChecked(visible)
        return layer

    def paletted(layer, entries):
        classes = [QgsPalettedRasterRenderer.Class(v, QColor(c), lbl)
                   for v, c, lbl in entries]
        layer.setRenderer(QgsPalettedRasterRenderer(
            layer.dataProvider(), 1, classes))

    def pseudocolor(layer, stops, vmin, vmax):
        shader = QgsRasterShader()
        ramp = QgsColorRampShader(vmin, vmax)
        ramp.setColorRampType(QgsColorRampShader.Interpolated)
        ramp.setColorRampItemList(
            [QgsColorRampShader.ColorRampItem(v, QColor(c), lbl)
             for v, c, lbl in stops])
        shader.setRasterShaderFunction(ramp)
        layer.setRenderer(QgsSingleBandPseudoColorRenderer(
            layer.dataProvider(), 1, shader))

    # ---------------------------------------------------------- basemap
    g_base = root.addGroup("Base")
    osm = QgsRasterLayer(
        "type=xyz&url=https://tile.openstreetmap.org/%7Bz%7D/%7Bx%7D/%7By%7D.png"
        "&zmax=19&zmin=0&http-header:referer=", "OpenStreetMap", "wms")
    if osm.isValid():
        add(osm, g_base)
    else:
        print("  ! OSM basemap layer invalid")

    for nm, f in (("Hillshade 90 m", "SriLanka_hillshade_90m.tif"),
                  ("Hillshade 30 m", "SriPada_hillshade_30m.tif")):
        l = QgsRasterLayer(os.path.join(DEMDIR, f), nm)
        if l.isValid():
            add(l, g_base, visible=False)

    for nm, f in (("DEM 90 m island", "SriLanka_DEM_90m.tif"),
                  ("DEM 30 m near field", "SriPada_DEM_30m.tif")):
        l = QgsRasterLayer(os.path.join(DEMDIR, f), nm)
        if l.isValid():
            add(l, g_base, visible=False)

    # ------------------------------------------------------- visibility
    g_vis = root.addGroup("Visibility")
    for nm, f, vis in (
            ("Observer height needed - island 90 m",
             "SriLanka_90m_class.tif", True),
            ("Observer height needed - near field 30 m",
             "SriPada_30m_class.tif", False)):
        l = QgsRasterLayer(os.path.join(VS, f), nm)
        if l.isValid():
            paletted(l, CLASS_COLORS)
            l.renderer().setOpacity(0.75)
            add(l, g_vis, visible=vis)

    for nm, f in (("How strongly it reads - island 90 m",
                   "SriLanka_90m_grade.tif"),
                  ("How strongly it reads - near field 30 m",
                   "SriPada_30m_grade.tif")):
        l = QgsRasterLayer(os.path.join(VS, f), nm)
        if l.isValid():
            paletted(l, GRADE_COLORS)
            l.renderer().setOpacity(0.8)
            add(l, g_vis, visible=False)

    for nm, f in (("Required height (m, continuous) 90 m",
                   "SriLanka_90m_required_height.tif"),
                  ("Required height (m, continuous) 30 m",
                   "SriPada_30m_required_height.tif")):
        l = QgsRasterLayer(os.path.join(VS, f), nm)
        if l.isValid():
            pseudocolor(l, [(-500, "#00441b", "500 m of headroom"),
                            (0, "#a1d99b", "just visible"),
                            (10, "#fed976", "10 m"),
                            (50, "#fd8d3c", "50 m"),
                            (150, "#bd0026", "150 m"),
                            (400, "#3f007d", "400 m+")], -500, 400)
            add(l, g_vis, visible=False)

    for nm, f in (("Apparent angular height (deg) 90 m",
                   "SriLanka_90m_angular_height_deg.tif"),
                  ("Distance to summit (km) 90 m",
                   "SriLanka_90m_distance_km.tif")):
        l = QgsRasterLayer(os.path.join(VS, f), nm)
        if l.isValid():
            add(l, g_vis, visible=False)

    # ------------------------------------------------------------ solar
    g_sol = root.addGroup("Sunrise geometry")
    for nm, f in (("Sun rises behind the peak - 1st date (day of year)",
                   "SriPada_30m_align_doy1.tif"),
                  ("Sun rises behind the peak - 2nd date (day of year)",
                   "SriPada_30m_align_doy2.tif")):
        p = os.path.join(SOL, f)
        if os.path.exists(p):
            l = QgsRasterLayer(p, nm)
            if l.isValid():
                pseudocolor(l, [(1, "#4575b4", "1 Jan"),
                                (80, "#91bfdb", "21 Mar"),
                                (172, "#fee090", "21 Jun"),
                                (266, "#fc8d59", "23 Sep"),
                                (355, "#d73027", "21 Dec")], 1, 365)
                add(l, g_sol, visible=False)

    for f in (sorted(os.listdir(SOL)) if os.path.isdir(SOL) else []):
        if f.startswith("SriPada_shadow_") and f.endswith(".tif"):
            l = QgsRasterLayer(os.path.join(SOL, f), f[:-4].replace("_", " "))
            if l.isValid():
                paletted(l, [(1, "#54278f", "terrain in shadow"),
                             (3, "#fdae6b", "shadow of the summit cone")])
                add(l, g_sol, visible=False)

    # ---------------------------------------------------------- vectors
    g_pt = root.addGroup("Points")
    sp = summit_geojson()
    l = QgsVectorLayer(sp, "Sri Pada summit", "ogr")
    if l.isValid():
        sym = QgsMarkerSymbol.createSimple(
            {"name": "triangle", "color": "#bd0026", "size": "6",
             "outline_color": "#ffffff", "outline_width": "0.6"})
        l.setRenderer(QgsSingleSymbolRenderer(sym))
        add(l, g_pt)

    pv = os.path.join(RES, "SriPada_places_visibility.geojson")
    if os.path.exists(pv):
        l = QgsVectorLayer(pv, "Places (with visibility attributes)", "ogr")
        if l.isValid():
            add(l, g_pt, visible=False)

    for nm, f in (("OSM places", "SriLanka_places.geojson"),
                  ("OSM buildings with height", "SriLanka_buildings.geojson"),
                  ("OSM peaks", "SriLanka_peaks.geojson")):
        p = os.path.join(OSMD, f)
        if os.path.exists(p):
            l = QgsVectorLayer(p, nm, "ogr")
            if l.isValid():
                add(l, g_pt, visible=False)

    out = os.path.join(ROOT, "Sri Pada Visibility.qgz")
    prj.write(out)
    print("wrote", out)
    qgs.exitQgis()


if __name__ == "__main__":
    print("=== hillshades ===")
    hillshades()
    print("=== qgis project ===")
    build_project()
