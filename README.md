# Sri Pada visibility study

Where in Sri Lanka can you see Sri Pada (Adam's Peak, 6.8094 N 80.4994 E), how
high would you have to be where you cannot, how large does it look, and where
does the rising sun come up from behind it.

Everything runs at **30 m**, COP30's native posting. The earlier two-tier
arrangement (30 m near field / 90 m island) is retired.

Run order under the QGIS python
(`"C:\Program Files\QGIS 3.40.8\bin\python-qgis-ltr.bat"`):

| step | script | what it does |
|---|---|---|
| 1 | `01_dem_prep.py` | pulls the 12 COP30 tiles, builds the mosaic VRT |
| 2 | `10_dem30.py` | island-wide 30 m DEM, land mask, hillshades |
| 3 | `11_visibility30.py` | the line-of-sight engine |
| 4 | `03_osm.py` | OSM places, buildings with height, peaks |
| 5 | `04_solar.py` | sunrise alignment + the peak's shadow |
| 6 | `13_analyse30.py` | land-masked statistics, per-town and per-building results |
| 7 | `12_web30.py` | the interactive map |
| 8 | `06_qgis.py` / `09_printmap.py` | QGIS project and the A2 sheet |
| 9 | `08_summary.py` | the numbers worth quoting |

## Deliverables

| what | where |
|---|---|
| **Public web app** (open in a browser, works offline, works on a phone) | `Web/SriPada_Visibility_Map.html` |
| QGIS project | `Sri Pada Visibility.qgz` |
| A2 print sheet | `SriPada_Visibility_Map_A2.pdf` / `.png` |
| 30 m rasters | `Viewshed/`, `Solar/` |
| 4,053 towns with distance, required height, sunrise dates | `Results/SriPada_places_visibility.csv` |
| Buildings whose roof wins the view | `Results/SriPada_buildings_visibility.csv` |

## The web app

Written for someone who has never heard of a viewshed. It answers one question —
*can you see Adam's Peak from here?* — in one of three ways:

- **Yes, from the ground.**
- **From about N metres up**, translated into storeys, so a roof or an upper
  floor is an answer people can act on.
- **No**, with the reason: either the specific height you would need, or that
  the curve of the earth hides it entirely.

Where the peak is visible it also gives the two mornings a year the sun rises
from directly behind the summit — the result that makes this more than a map.

Other things built for public use: search across all 4,004 named settlements
with the answer colour-coded before you click; a locate-me button; a shareable
link that restores the exact spot, bearing and eye height; a first-run
explainer; and a mobile layout with a bottom sheet. Deliberately shown in plain
language — "how big it looks", not "apparent angular height"; "earth curvature
hides 351 m of it", not a refraction coefficient.

The second tab draws the skyline you would actually see toward the peak, with
a sunrise-light mode.

## Method

**Line of sight.** Rays are cast *from the summit* on 24,000 azimuths, sampling
the DEM every 30 m out to 190 km. Along each ray the running maximum of the
vertical angle subtended at the summit,

    alpha(r) = ( z(r) - c(r) - Zs ) / r

is the grazing angle of the highest obstruction so far. A cell at horizontal
distance `d` then needs its surface to reach `z_req = Zs + A* d + c(d)`.
Because line of sight is reciprocal, one outward sweep answers every observer
at once instead of one viewshed per observer.

**Earth curvature and atmospheric refraction are both modelled**, everywhere:

    c(d) = (1 - k) d^2 / (2R),   k = 0.13,   R = 6,371,008.8 m

k is the standard terrestrial refraction coefficient — light bends toward the
earth, so effective curvature is only (1−k) of geometric. It is applied to
every terrain sample along every ray *and* to the target cell. This is not a
cosmetic correction: it is a 351 m drop at Colombo's 72 km and 2,731 m at
200 km, and it moves the sea-level horizon for a 2,192 m peak from 172 km
(k = 0) to 184 km. Rays stop at 190 km because past that the curvature drop
alone exceeds the summit's height, so nothing at ground level can see it.

**The master layer is signed.** `required_height = z_req − z_ground` is
positive where it is the metres of observer height you need (the tall-building
question) and negative where it is the metres of headroom you have. Both halves
of the question are the same number.

**Grade.** Apparent angular height of the summit, binned: prominent ≥ 2°,
clear 1–2°, distinct 0.5–1°, faint 0.2–0.5°, marginal < 0.2°. This is what
separates a real sighting from bare geometry.

**Sunrise alignment.** The summit sits at a fixed point in an observer's sky
(bearing A, apparent altitude θ). The sun crosses it when

    sin(dec) = sin(lat) sin(h) + cos(lat) cos(h) cos(A)

with h = θ minus astronomical refraction. Declination sweeps ±23.44° once a
year, so a solvable cell gets exactly two dates. Since
`sin(H) = −cos(h) sin(A) / cos(dec)`, only observers **west** of the peak meet
the sun on its rising branch — that one inequality carves the corridor out of
the map. Recovered bearing window 66.4°–114.9° against a theoretical
66.4°–113.6°; the slight overshoot is real, high observers see the summit at
positive altitude.

**Shadow.** Terrain is marched along the solar bearing keeping the running
maximum of `z − s·tan(alt)`; the index of the *blocking* cell is carried
through, so the summit cone's own shadow is separated from every other ridge's.

## Headline results (30 m, land only, 65,896 km²)

| observer height | new area | cumulative | share of land |
|---|---|---|---|
| ground level | 7,965 km² | 7,965 | 12.1% |
| ≤ 10 m (3 storeys) | 4,667 | 12,632 | 19.2% |
| ≤ 30 m (10 storeys) | 3,230 | 15,862 | 24.1% |
| ≤ 60 m (20 storeys) | 2,555 | 18,418 | 27.9% |
| ≤ 100 m (33 storeys) | 2,100 | 20,517 | 31.1% |
| ≤ 150 m (50 storeys) | 1,794 | 22,311 | 33.9% |
| never, at any height | | 43,584 | 66.1% |

Furthest ground-level sighting **by how strongly the peak reads** — the honest
framing, because bare geometry is misleading at range:

| grade | furthest place | km | apparent height |
|---|---|---|---|
| prominent ≥ 2° | Gelanigama | 56.0 | 2.01° |
| clear 1–2° | Unawatuna | 91.6 | 1.01° |
| distinct 0.5–1° | Battulu Oya | 125.4 | 0.51° |
| faint 0.2–0.5° | Kumbukwewa | 151.4 | 0.21° |
| marginal < 0.2° | Kalpitiya | 177.4 | 0.01° |

Kalpitiya's 177 km is real geometry but not a sighting: 0.01° is a fiftieth of
the moon's width. Use the grade layer, never the binary viewshed, for anything
about what people actually see.

Of 4,053 towns, **705 see it from the ground** and 1,063 more would from a
rooftop. Of 180,182 OSM buildings carrying a height or level tag, **11,266
win the view from the roof when their footprint does not** — Colombo sits right
on the threshold, which is why the Lotus Tower needs just 4.8 m at 72 km. The
rising sun comes up from directly behind the summit, seen from Colombo, on
**19 February and 22 October**; from Kalutara on 6 June and 6 July, straddling
the solstice as its 67° bearing demands.

## Caveats that matter

- **The DEM is a DSM.** Canopy is included: it blocks realistically, but it
  also lifts forested observers to treetop height, so ground-level visibility
  in forest is optimistic. FABDEM is the forest-removed alternative.
- **COP30 truncates the summit to 2,192 m** against a surveyed 2,243 m — a
  30 m posting cannot hold a sharp cone. Summit elevation is therefore set
  explicitly rather than read off the grid. Results are mildly conservative.
- **Land mask** is `z != 0` (Copernicus codes sea as exact zero), recovering
  65,896 km² against a surveyed 65,610 km², 0.4% out. All quoted areas are land
  only — the raw analysis box is mostly ocean.
- **Haze is not modelled.** Typical tropical visual range is 20–60 km.
- The web map's click readout uses a ~220 m resample for terrain elevation and
  ridge shape; the visibility verdict always comes from the 30 m raster.
- Leaflet is inlined in the HTML, so the map needs no network at all. The
  default basemap is baked from this study's own hillshade; the online
  providers (Carto, Esri, OpenTopoMap, OSM) are optional extras.

## See also

`..\SriLankaPeaks\` — the same machinery generalised to all 1,237 peaks on the
island, with visibility solved in both directions.
