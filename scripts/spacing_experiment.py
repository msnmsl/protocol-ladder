# -*- coding: utf-8 -*-
"""Regenerate the sampling-density manipulation behind Table S4.

Part 1 is observational: 100 random draws of 40 compositions x 5 temperature
points, regressing the paired loss on realised within-composition spacing.
Part 2 is the controlled comparison: 30 independent draws per condition. The
two conditions apply different placement rules (consecutive vs quantile-spaced)
to independently drawn compositions, so the comparison is two-sample rather
than paired: the placement rule is the design variable, not the identity of
the compositions drawn.

Needs only the public gpvisc CSV. Set DSTAR_GPVISC_CSV to relocate it and
DSTAR_OUT to redirect the two output files.
"""
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import GroupKFold, KFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paths import GPVISC_CSV, OUT

warnings.filterwarnings("ignore")

SEED, NS, NG, PG = 42, 5, 40, 5
FEATS = ["sio2", "al2o3", "cao", "mgo", "R", "invT"]


def r2(y, p):
    return 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)


def cv(X, y, g=None):
    p = np.full(len(y), np.nan)
    it = (KFold(NS, shuffle=True, random_state=SEED).split(X) if g is None
          else GroupKFold(NS).split(X, y, groups=g))
    for tr, te in it:
        m = ExtraTreesRegressor(n_estimators=200, random_state=SEED, n_jobs=-1).fit(X[tr], y[tr])
        p[te] = m.predict(X[te])
    return r2(y, p)


def loss_of(sub):
    assert sub["Name"].nunique() >= NS, f"only {sub['Name'].nunique()} groups"
    a = cv(sub[FEATS].values, sub["lnEta"].values)
    b = cv(sub[FEATS].values, sub["lnEta"].values, sub["Name"].values)
    return (a - b) if (np.isfinite(a) and np.isfinite(b)) else np.nan


raw = pd.read_csv(GPVISC_CSV)
d = raw.dropna(subset=["viscosity", "T", "Name"]).copy()
d = d[np.isfinite(d["viscosity"]) & np.isfinite(d["T"])]
d = d[d["T"] >= 1000.0]
sz = d.groupby("Name").size()
d = d[d["Name"].isin(sz[sz >= 5].index)].copy()
d["lnEta"] = d["viscosity"] * math.log(10.0)
d = d[d["lnEta"] <= 12.0].reset_index(drop=True)
d["invT"] = 1000.0 / d["T"]
d["R"] = d["cao"] / d["sio2"].replace(0, np.nan)
d["R"] = d["R"].fillna(d["R"].median())
d = d.sort_values(["Name", "T"]).reset_index(drop=True)
print(f"base: {len(d)} points / {d['Name'].nunique()} groups", flush=True)

by = {nm: v for nm, v in d.groupby("Name")}   # keep the global index; do NOT reset_index
names5 = [k for k, v in by.items() if len(v) >= PG]
names8 = [k for k, v in by.items() if len(v) >= 8]
print(f"groups with >={PG} points: {len(names5)}; with >=8: {len(names8)}", flush=True)

rng = np.random.default_rng(SEED)
rows = []

print("\n[Part 1] random 40x5 draws, regress loss on realised spacing", flush=True)
rec = []
for it in range(100):
    pick = rng.choice(names5, NG, replace=False)
    idx, sp = [], []
    for nm in pick:
        v = by[nm]
        take = rng.choice(len(v), PG, replace=False)
        s = v.iloc[take].sort_values("T")
        idx.extend(s.index.tolist())
        tt = s["T"].values
        sp.append(float(np.mean(np.diff(tt))))
    sub = d.loc[idx]
    L = loss_of(sub)
    if np.isfinite(L):
        rec.append(dict(draw=it, spacing=float(np.mean(sp)), loss=L))
p1 = pd.DataFrame(rec)
r, pv = stats.spearmanr(p1["spacing"], p1["loss"])
sl, ic, rr, pp, se = stats.linregress(p1["spacing"], p1["loss"])
print(f"  n={len(p1)} draws  spacing {p1['spacing'].min():.0f}-{p1['spacing'].max():.0f} K "
      f"(median {p1['spacing'].median():.0f})", flush=True)
print(f"  Spearman rho = {r:+.3f} (p={pv:.4f}); slope = {sl:+.5f} per K (R2={rr**2:.3f})", flush=True)
rows.append(dict(part=1, condition="all_draws", n=len(p1),
                 spacing=float(p1["spacing"].median()), loss=float(p1["loss"].mean()),
                 sd=float(p1["loss"].std()), rho=float(r), p=float(pv)))

print("\n[Part 2] same sampling frame, tight vs wide 5-point placement", flush=True)
res = {}
for cond in ("tight", "wide"):
    out, sp = [], []
    for it in range(30):
        pick = rng.choice(names8, NG, replace=False)
        idx = []
        for nm in pick:
            v = by[nm]
            n = len(v)
            if cond == "tight":
                start = rng.integers(0, n - PG + 1)
                s = v.iloc[start:start + PG]
            else:
                pos = np.linspace(0, n - 1, PG).round().astype(int)
                s = v.iloc[pos]
            idx.extend(s.index.tolist())
            tt = np.sort(s["T"].values)
            sp.append(float(np.mean(np.diff(tt))))
        sub = d.loc[idx]
        L = loss_of(sub)
        if np.isfinite(L):
            out.append(L)
    res[cond] = (float(np.mean(out)), float(np.std(out)), float(np.mean(sp)))
    print(f"  {cond:5s}: realised spacing {np.mean(sp):5.1f} K  loss = {np.mean(out):+.4f} "
          f"(sd {np.std(out):.3f}, k={len(out)})", flush=True)
    rows.append(dict(part=2, condition=cond, n=len(out), spacing=float(np.mean(sp)),
                     loss=float(np.mean(out)), sd=float(np.std(out)), rho=np.nan, p=np.nan))

if "tight" in res and "wide" in res:
    diff = res["tight"][0] - res["wide"][0]
    se_d = math.sqrt(res["tight"][1] ** 2 / 30 + res["wide"][1] ** 2 / 30)
    print(f"\n  tight - wide = {diff:+.4f}  (SE {se_d:.4f}, t = {diff/se_d:+.2f})", flush=True)

out_df = pd.DataFrame(rows)
out_df.to_csv(OUT / "spacing_experiment.csv", index=False, encoding="utf-8-sig")
p1.to_csv(OUT / "spacing_part1_draws.csv", index=False, encoding="utf-8-sig")
print(f"\nwrote {OUT / 'spacing_experiment.csv'}")
print(out_df.round(4).to_string(index=False))
