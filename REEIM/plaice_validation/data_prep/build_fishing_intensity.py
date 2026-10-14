"""Relative fishing intensity for North Sea plaice, 1850-2020, scaled to a maximum of 1.

1850-1892: assumption. Linear ramp from 0 in RAMP_START to the 1893 value (sail trawling,
           then steam trawling from the 1880s; no F estimates before 1893).

1893-1956: Rijnsdorp & Millner (1996) Table 1 F from c.p.u.e., rescaled to the ICES 2024
           F(2-6) scale with the 1957-66 overlap (factor 0.5116; see
           data/plaice_F_pre1957_splice_derived.csv).
War years: the source has no F for 1915-18 and 1939-45. Assumption: effort fell to
           WAR_FRACTION of the mean F of the four pre-war years (British steam-trawler effort
           fell 97% from 1938 to 1941; no WWI figure found). 1939 is the mean of 1938 and 1940.
1957-2020: ICES SAG 2024 assessment F (ages 2-6).
Usage (from the repo root):
  python3 REEIM/plaice_validation/data_prep/build_fishing_intensity.py REEIM/plaice_validation/data_prep/sources REEIM/input/plaice/fishing_intensity.csv [war_fraction=0.1] [terminal_year=2020]
"""
import sys, pandas as pd, numpy as np
d, out = sys.argv[1], sys.argv[2]
WAR_FRACTION = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
TERM = int(sys.argv[4]) if len(sys.argv) > 4 else 2020

old = pd.read_csv(f"{d}/plaice_F_pre1957_splice_derived.csv", comment="#")
ices = pd.read_csv(f"{d}/ices_sag_ple.27.420_2024.csv", comment="#")
F = dict(zip(old.year, old.F_ICESscale_cpue))
F.update(dict(zip(ices.Year, ices.F)))
for war, pre in [(range(1915, 1919), range(1910, 1914)), (range(1940, 1946), range(1935, 1939))]:
    lvl = np.mean([F[y] for y in pre]) * WAR_FRACTION
    for y in war: F[y] = lvl
F[1939] = (F[1938] + F[1940]) / 2
RAMP_START = 1850
for y in range(RAMP_START, 1893):
    F[y] = F[1893] * (y - RAMP_START) / (1893 - RAMP_START)
years = list(range(RAMP_START, TERM + 1))
f = np.array([F[y] for y in years], float)
assert not np.isnan(f).any(), [y for y, v in zip(years, f) if np.isnan(v)]
pd.DataFrame({"year": years, "intensity": np.round(f / f.max(), 3)}).to_csv(out, index=False)
print(f"max F {f.max():.3f} in {years[f.argmax()]}; terminal {TERM} relative = {f[-1]/f.max():.3f}")
