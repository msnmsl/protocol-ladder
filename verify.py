# -*- coding: utf-8 -*-
"""Check every number the protocol paper quotes against the shipped CSVs.

Exit code is non-zero if any check fails, so this can gate a release.

    python verify.py
"""
from __future__ import annotations

import math
import sys

import pandas as pd
from scipy import stats

from paths import OUT  # honours DSTAR_OUT

fails: list[str] = []
checked = 0


def check(label: str, got, want, tol=5e-4) -> None:
    global checked
    checked += 1
    ok = abs(float(got) - float(want)) <= tol if isinstance(want, (int, float)) else got == want
    print(f"  {'ok  ' if ok else 'FAIL'} {label:52s} {got!r:>10}  (expected {want})")
    if not ok:
        fails.append(label)


print("=== plant record, 22 configurations x 4 tiers ===")
d = pd.read_csv(OUT / "g4_plant_ladder.csv", encoding="utf-8-sig")
check("rows", len(d), 88)
check("configurations", d.groupby(["模型", "方案"]).ngroups, 22)
med = d.pivot_table(index=["模型", "方案"], columns="协议", values="R2").median()
for t, want in [("R1", 0.9773), ("R2", 0.7274), ("R3", 0.0534), ("R4", -1.5302)]:
    check(f"median R2 {t}", round(float(med[t]), 4), want)
et = d[(d["模型"] == "ExtraTrees") & (d["方案"] == "F3摩尔+碱度")].set_index("协议")
for t, want in [("R1", 0.0000), ("R2", 0.2163), ("R3", 0.4804), ("R4", 0.5575)]:
    check(f"ET-F3 RMSE {t}", round(float(et.loc[t, "rmse_ln"]), 4), want)

print("\n=== wide-range database ===")
g = pd.read_csv(OUT / "g4_gpvisc_ladder.csv", encoding="utf-8-sig")
check("rows", len(g), 88)
ge = g[(g["模型"] == "ExtraTrees") & (g["方案"] == "F3-mol+R+T")].set_index("协议")
for t, want in [("R1", 0.4312), ("R2", 0.5381), ("R3", 0.6262), ("R4", 0.9255)]:
    check(f"ET-F3 RMSE {t}", round(float(ge.loc[t, "rmse_ln"]), 4), want)

print("\n=== literature compilation ===")
c = pd.read_csv(OUT / "g4_compilation_ladder.csv", encoding="utf-8-sig")
check("rows", len(c), 88)
ce = c[(c["模型"] == "ExtraTrees") & (c["方案"] == "F1-4氧化物质分数")].set_index("协议")
check("ET-F1 R2 T2", round(float(ce.loc["R2", "R2"]), 3), 0.631)
check("ET-F1 R2 T3", round(float(ce.loc["R3", "R2"]), 3), 0.148)

print("\n=== paired loss T2 -> T3 (median over 22 configurations) ===")
for name, frame in [("plant", d), ("wide-range", g), ("compilation", c)]:
    piv = frame.pivot_table(index=["模型", "方案"], columns="协议", values="R2")
    loss = (piv["R2"] - piv["R3"]).median()
    label = {"plant": 0.7644, "wide-range": 0.0086, "compilation": 0.7724}[name]
    check(f"{name} median paired loss", round(float(loss), 4), label)

print("\n=== sampling-density manipulation (Table S4) ===")
s4 = pd.read_csv(OUT / "spacing_experiment.csv", encoding="utf-8-sig")
tv = s4[s4["condition"] == "tight"].iloc[0]
wv = s4[s4["condition"] == "wide"].iloc[0]
T = 1e-4  # tighter than the default for four-decimal values
check("tight realised spacing (K)", round(float(tv["spacing"]), 1), 32.0, T)
check("wide realised spacing (K)", round(float(wv["spacing"]), 1), 79.5, T)
check("tight repeats", int(tv["n"]), 30, T)
check("wide repeats", int(wv["n"]), 30, T)
check("tight paired loss", round(float(tv["loss"]), 4), 0.2100, T)
check("wide paired loss", round(float(wv["loss"]), 4), 0.1185, T)
check("tight sd", round(float(tv["sd"]), 4), 0.0849, T)
check("wide sd", round(float(wv["sd"]), 4), 0.0435, T)
d4 = float(tv["loss"]) - float(wv["loss"])
check("tight - wide difference", round(d4, 4), 0.0915, T)
se4 = math.sqrt(float(tv["sd"]) ** 2 / float(tv["n"]) + float(wv["sd"]) ** 2 / float(wv["n"]))
check("SE of difference", round(se4, 4), 0.0174, T)
check("t statistic", round(d4 / se4, 2), 5.25, T)
p4 = pd.read_csv(OUT / "spacing_part1_draws.csv", encoding="utf-8-sig")
check("Part 1 draws", len(p4), 100, T)
rho4, pv4 = stats.spearmanr(p4["spacing"], p4["loss"])
check("Part 1 Spearman rho", round(float(rho4), 3), -0.304, T)
check("Part 1 Spearman p", round(float(pv4), 4), 0.0021, T)
sl4, ic4, rr4, pp4, se_r4 = stats.linregress(p4["spacing"], p4["loss"])
check("Part 1 slope per K", round(float(sl4), 5), -0.00226, 1e-5)
check("Part 1 R2", round(float(rr4) ** 2, 3), 0.109, T)
check("Part 1 spacing min (K)", round(float(p4["spacing"].min()), 1), 56.9, T)
check("Part 1 spacing max (K)", round(float(p4["spacing"].max()), 1), 91.1, T)
check("Part 1 spacing median (K)", round(float(p4["spacing"].median()), 1), 74.4, T)
check("Part 1 intercept", round(float(ic4), 3), 0.293, T)
pred9 = float(ic4) + float(sl4) * 9.0
check("extrapolated loss at 9 K", round(pred9, 3), 0.273, T)
check("share of plant loss 0.7644", round(pred9 / 0.7644, 2), 0.36, T)

print(f"\n{checked - len(fails)}/{checked} checks passed")
if fails:
    print("FAILED:", *fails, sep="\n  ")
    sys.exit(1)
print("every quoted number traces to a shipped file")
