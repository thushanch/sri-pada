"""
Sri Pada visibility study - Step 3: OSM layers.

Pulls three things from Overpass, island-wide:

  places     cities / towns / villages, used as the reporting points and for
             labelling the map
  buildings  anything carrying height or building:levels, so real rooftops can
             be tested against the required-height surface
  peaks      natural=peak, for cross-checking the DEM and for context

Note on coverage: height tagging in Sri Lanka is sparse. The building layer
validates the required-height raster at real high-rises, it does not map the
country's building stock.
"""
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request

ROOT = r"C:\Users\thush\OneDrive\Desktop\SriPada"
OUT = os.path.join(ROOT, "OSM")

ENDPOINTS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
    "https://overpass.private.coffee/api/interpreter",
]

# Sri Lanka bounding box (S, W, N, E)
BBOX = "5.7,79.4,10.0,82.0"

QUERIES = {
    "places": f"""
[out:json][timeout:300];
(
  node["place"~"^(city|town|village|suburb)$"]({BBOX});
);
out body;
""",
    "buildings": f"""
[out:json][timeout:300];
(
  way["building"]["height"]({BBOX});
  way["building"]["building:levels"]({BBOX});
  relation["building"]["height"]({BBOX});
  relation["building"]["building:levels"]({BBOX});
);
out center tags;
""",
    "peaks": f"""
[out:json][timeout:300];
(
  node["natural"="peak"]({BBOX});
);
out body;
""",
}


def overpass(q, tries=3):
    data = urllib.parse.urlencode({"data": q}).encode()
    last = None
    for attempt in range(tries):
        for ep in ENDPOINTS:
            try:
                req = urllib.request.Request(
                    ep, data=data,
                    headers={"User-Agent": "SriPada-visibility-study/1.0"})
                with urllib.request.urlopen(req, timeout=400) as r:
                    return json.loads(r.read().decode())
            except Exception as e:                      # noqa: BLE001
                last = f"{ep}: {e}"
                print("   retry:", last, flush=True)
                time.sleep(5)
    raise RuntimeError(f"overpass failed: {last}")


def to_geojson(elements, kind):
    feats = []
    for el in elements:
        if el["type"] == "node":
            lon, lat = el.get("lon"), el.get("lat")
        else:
            c = el.get("center") or {}
            lon, lat = c.get("lon"), c.get("lat")
        if lon is None or lat is None:
            continue
        tags = el.get("tags", {}) or {}
        props = {"osm_id": el["id"], "osm_type": el["type"]}
        for k in ("name", "name:en", "place", "natural", "ele", "height",
                  "building:levels", "building", "population"):
            if k in tags:
                props[k.replace(":", "_")] = tags[k]

        if kind == "buildings":
            h = None
            raw = tags.get("height")
            if raw:
                try:
                    h = float(str(raw).lower().replace("m", "").strip())
                except ValueError:
                    h = None
            if h is None and tags.get("building:levels"):
                try:
                    h = float(str(tags["building:levels"]).split(";")[0]) * 3.2
                except ValueError:
                    h = None
            if h is None:
                continue
            props["height_m"] = round(h, 1)
            props["height_src"] = "height" if tags.get("height") else "levels x3.2"

        feats.append({"type": "Feature",
                      "geometry": {"type": "Point", "coordinates": [lon, lat]},
                      "properties": props})
    return {"type": "FeatureCollection", "features": feats}


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for kind, q in QUERIES.items():
        print(f"querying {kind} ...", flush=True)
        res = overpass(q)
        gj = to_geojson(res.get("elements", []), kind)
        p = os.path.join(OUT, f"SriLanka_{kind}.geojson")
        with open(p, "w", encoding="utf-8") as f:
            json.dump(gj, f, ensure_ascii=False)
        n = len(gj["features"])
        print(f"  {kind}: {n} features -> {os.path.basename(p)}", flush=True)
        if kind == "buildings" and n:
            hs = sorted((f["properties"]["height_m"] for f in gj["features"]),
                        reverse=True)
            print(f"  tallest tagged: {hs[:10]}")
        time.sleep(3)
