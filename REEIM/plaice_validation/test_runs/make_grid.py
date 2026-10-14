"""Write SLiM commands for a plaice grid over S_LT x w at a fixed H (other parameters at parameters.csv defaults).

Usage (from the repository root, or a copy of it holding REEIM/REEIM.slim and REEIM/input/):
  python3 REEIM/plaice_validation/test_runs/make_grid.py <out_dir> [H=0.078] [reps=2] [Tdelta=10000] [set=plaice]
Writes <out_dir>/cmds.txt and <out_dir>/design.csv in the make_test_runs.py format, so
analyse_test_runs.R summarises it the same way.
"""
import csv, os, sys
out_dir = os.path.abspath(sys.argv[1]); H = float(sys.argv[2]) if len(sys.argv) > 2 else 0.078
reps = int(sys.argv[3]) if len(sys.argv) > 3 else 2; Tdelta = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
st = sys.argv[5] if len(sys.argv) > 5 else "plaice"
root = os.getcwd(); model = os.path.join(root, "REEIM", "REEIM.slim"); inp = os.path.join(root, "REEIM", "input", st)
pr = {r["parameter"]: (float(r["lower"]), float(r["upper"])) for r in csv.DictReader(open(os.path.join(inp, "priors.csv")))}
pa = {r[0]: r[1] for r in csv.reader(open(os.path.join(inp, "parameters.csv"))) if r and not r[0].startswith("#")}
lv = lambda k: [pr[k][0] + i * (pr[k][1] - pr[k][0]) / 4 for i in range(5)]
fol = os.path.join(out_dir, "out", st); os.makedirs(fol, exist_ok=True)
f, Shmax = float(pa["f"]), float(pa["Shmax"]); i = 0
with open(os.path.join(out_dir, "cmds.txt"), "w") as out, open(os.path.join(out_dir, "design.csv"), "w") as d:
    d.write("set,rep,H,f,w,Shmax,SLT\n")
    for SLT in lv("SLT"):
        for w in lv("w"):
            for r in range(reps):
                i += 1
                args = f"-d H={H} -d f={f} -d w={round(w, 3)} -d Shmax={round(Shmax)} -d SLT={round(SLT)}"
                out.write(f"slim -s {2000 + i} -d \"folder='{fol}'\" -d \"inputDir='{inp}'\" -d Tdelta={Tdelta} -d rep={i} {args} {model} > /dev/null 2>&1\n")
                d.write(f"{st},{i},{H},{f},{round(w, 3)},{round(Shmax)},{round(SLT)}\n")
