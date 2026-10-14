"""Write SLiM commands for plaice test runs: Latin-hypercube draws from each input set's priors.

Usage (from the repository root, or a copy of it holding REEIM/REEIM.slim and REEIM/input/):
  python3 REEIM/plaice_validation/test_runs/make_test_runs.py <out_dir> [n_per_set=24] [Tdelta=10000]
Writes <out_dir>/cmds.txt (one SLiM call per line) and <out_dir>/design.csv. Run e.g. with
  cat <out_dir>/cmds.txt | parallel -j 9      (or: xargs -P 9 -n 1 bash over one script per line)
Outputs go to <out_dir>/out/<set>/. Summarise with test_runs/analyse_test_runs.R.
"""
import csv, os, random, sys
out_dir = os.path.abspath(sys.argv[1]); n = int(sys.argv[2]) if len(sys.argv) > 2 else 24
Tdelta = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
random.seed(1)
root = os.getcwd(); model = os.path.join(root, "REEIM", "REEIM.slim")
with open(os.path.join(out_dir, "cmds.txt"), "w") as out, open(os.path.join(out_dir, "design.csv"), "w") as d:
    d.write("set,rep,H,f,w,Shmax,SLT\n")
    for st in ["plaice", "plaice_gfree"]:
        inp = os.path.join(root, "REEIM", "input", st)
        pr = {r["parameter"]: (float(r["lower"]), float(r["upper"])) for r in csv.DictReader(open(os.path.join(inp, "priors.csv")))}
        cols = {k: [(i + random.random()) / n for i in range(n)] for k in pr}
        for k in cols: random.shuffle(cols[k])
        fol = os.path.join(out_dir, "out", st); os.makedirs(fol, exist_ok=True)
        for i in range(n):
            th = {k: pr[k][0] + cols[k][i] * (pr[k][1] - pr[k][0]) for k in pr}
            th = {"H": round(th["H"], 3), "f": round(th["f"], 4), "w": round(th["w"], 3), "Shmax": round(th["Shmax"]), "SLT": round(th["SLT"])}
            args = " ".join(f"-d {k}={v}" for k, v in th.items())
            out.write(f"slim -s {1000 + i} -d \"folder='{fol}'\" -d \"inputDir='{inp}'\" -d Tdelta={Tdelta} -d rep={i + 1} {args} {model} > /dev/null 2>&1\n")
            d.write(f"{st},{i + 1},{th['H']},{th['f']},{th['w']},{th['Shmax']},{th['SLT']}\n")
