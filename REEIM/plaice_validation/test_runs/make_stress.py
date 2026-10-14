"""Write SLiM commands for the plaice maximum-response (stress) test: how much L_T decline can the model produce at all?

Usage (from the repository root, or a copy of it holding REEIM/REEIM.slim and REEIM/input/):
  python3 REEIM/plaice_validation/test_runs/make_stress.py <out_dir> [reps=2] [Tdelta=10000] [set=plaice]
Design: H {0.078, 0.2} x SLT {324, 420, 500} x w {0.02, 0.1} x genetics {default, more, extreme} x reps; f = 0.1, Shmax = 117.
Genetics levels override the global sd, u, Ve (all three traits): default = parameters.csv; more = Ve 0.02, sd 0.1;
extreme = Ve 0.01, sd 0.1, u 1e-5. Writes <out_dir>/cmds.txt and <out_dir>/design.csv (make_test_runs.py format plus a
'gen' column) so analyse_test_runs.R summarises it the same way.
"""
import csv, os, sys
out_dir = os.path.abspath(sys.argv[1]); reps = int(sys.argv[2]) if len(sys.argv) > 2 else 2
Tdelta = int(sys.argv[3]) if len(sys.argv) > 3 else 10000; st = sys.argv[4] if len(sys.argv) > 4 else "plaice"
root = os.getcwd(); model = os.path.join(root, "REEIM", "REEIM.slim"); inp = os.path.join(root, "REEIM", "input", st)
pa = {r[0]: r[1] for r in csv.reader(open(os.path.join(inp, "parameters.csv"))) if r and not r[0].startswith("#")}
gen = {"default": {"Ve": float(pa["Ve"]), "sd": float(pa["sd"]), "u": float(pa["u"])},
       "more": {"Ve": 0.02, "sd": 0.1, "u": float(pa["u"])},
       "extreme": {"Ve": 0.01, "sd": 0.1, "u": 1e-5}}
fol = os.path.join(out_dir, "out", st); os.makedirs(fol, exist_ok=True)
f, Shmax = 0.1, round(float(pa["Shmax"])); i = 0
with open(os.path.join(out_dir, "cmds.txt"), "w") as out, open(os.path.join(out_dir, "design.csv"), "w") as d:
    d.write("set,rep,H,f,w,Shmax,SLT,gen,Ve,sd,u\n")
    for H in (0.078, 0.2):
        for SLT in (324, 420, 500):
            for w in (0.02, 0.1):
                for gk, g in gen.items():
                    for r in range(reps):
                        i += 1
                        args = (f"-d H={H} -d f={f} -d w={w} -d Shmax={Shmax} -d SLT={SLT} "
                                f"-d Ve={g['Ve']} -d sd={g['sd']} -d u={g['u']}")
                        out.write(f"slim -s {3000 + i} -d \"folder='{fol}'\" -d \"inputDir='{inp}'\" -d Tdelta={Tdelta} -d rep={i} {args} {model} > /dev/null 2>&1\n")
                        d.write(f"{st},{i},{H},{f},{w},{Shmax},{SLT},{gk},{g['Ve']},{g['sd']},{g['u']}\n")
