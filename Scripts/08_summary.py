"""
Sri Pada visibility study - Step 8: the numbers worth quoting.

Separates geometric visibility from perceptible visibility. A line of sight
that clears the terrain by a metre at 177 km is real geometry but not a real
sighting: at that range the summit stands under a hundredth of a degree above
the horizon, roughly a fiftieth of the moon's width, and haze closes it out
long before the eye does.
"""
import csv
import os

RES = r"C:\Users\thush\OneDrive\Desktop\SriPada\Results"

GRADE = [(2.0, "prominent"), (1.0, "clear"), (0.5, "distinct"),
         (0.2, "faint"), (0.0, "marginal")]

CITIES = ["Colombo", "Dehiwala-Mount Lavinia", "Sri Jayawardenepura Kotte",
          "Negombo", "Kandy", "Galle", "Matara", "Ratnapura", "Nuwara Eliya",
          "Badulla", "Kurunegala", "Anuradhapura", "Polonnaruwa",
          "Trincomalee", "Batticaloa", "Jaffna", "Hambantota", "Puttalam",
          "Kalpitiya", "Chilaw", "Kalutara", "Avissawella", "Hatton",
          "Bandarawela", "Ella", "Matale", "Gampaha", "Vavuniya"]


def load():
    with open(os.path.join(RES, "SriPada_places_visibility.csv"),
              encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def fnum(r, k):
    try:
        return float(r[k])
    except (ValueError, KeyError, TypeError):
        return None


if __name__ == "__main__":
    rows = load()
    ground = [r for r in rows if r["visibility"] == "ground level"]

    print("=" * 74)
    print("FURTHEST GROUND-LEVEL SIGHTING, BY HOW STRONGLY THE PEAK READS")
    print("=" * 74)
    print(f"{'threshold':<26}{'furthest place':<26}{'km':>7}{'deg':>8}")
    for thr, name in GRADE:
        sub = [r for r in ground
               if (fnum(r, "angular_height_deg") or 0) >= thr]
        if not sub:
            continue
        far = max(sub, key=lambda r: float(r["distance_km"]))
        lbl = f"{name} (>={thr}deg)"
        print(f"{lbl:<26}{(far['name'] or '(unnamed)')[:24]:<26}"
              f"{float(far['distance_km']):>7.1f}"
              f"{fnum(far, 'angular_height_deg'):>8.2f}")
    print(f"\n{'any (>0deg, geometric only)':<26}"
          f"{max(ground, key=lambda r: float(r['distance_km']))['name'][:24]:<26}"
          f"{max(float(r['distance_km']) for r in ground):>7.1f}")

    n_al = sum(1 for r in rows if r["sunrise_align_1"])
    print(f"\nplaces where the rising sun comes up from behind the peak: {n_al}"
          f"  of {len(rows)}")

    print("\n" + "=" * 74)
    print("MAJOR TOWNS")
    print("=" * 74)
    print(f"{'town':<28}{'km':>7}{'bearing':>9}{'needs':>9}{'apparent':>10}"
          f"  {'sun behind peak':<22}")
    seen = set()
    for c in CITIES:
        m = [r for r in rows if r["name"] == c]
        if not m:
            continue
        r = min(m, key=lambda x: float(x["distance_km"]))
        if r["name"] in seen:
            continue
        seen.add(r["name"])
        need = float(r["required_height_m"])
        needs = "ground" if need <= 1.7 else (
            f"{need:.0f} m" if need <= 300 else "-")
        ang = fnum(r, "angular_height_deg")
        al = (f"{r['sunrise_align_1']} & {r['sunrise_align_2']}"
              if r["sunrise_align_1"] else "-")
        print(f"{c[:27]:<28}{float(r['distance_km']):>7.1f}"
              f"{float(r['bearing_to_summit_deg']):>8.0f}d{needs:>9}"
              f"{(f'{ang:.2f}d' if ang else '-'):>10}  {al:<22}")

    print("\n" + "=" * 74)
    print("BUILDINGS THAT WIN THE VIEW FROM THE ROOF")
    print("=" * 74)
    with open(os.path.join(RES, "SriPada_buildings_visibility.csv"),
              encoding="utf-8-sig") as f:
        b = list(csv.DictReader(f))
    print(f"total: {len(b):,}")
    named = [r for r in b if r["name"]]
    print(f"named: {len(named):,}")
    by_need = sorted(b, key=lambda r: -float(r["required_height_m"]))
    print("\nthe hardest-won views (tallest requirement actually met):")
    for r in by_need[:10]:
        print(f"  {(r['name'] or '(unnamed)')[:34]:<36}"
              f"{float(r['building_height_m']):>6.0f} m roof   needs "
              f"{float(r['required_height_m']):>6.1f} m   @"
              f"{float(r['distance_km']):>6.1f} km")
