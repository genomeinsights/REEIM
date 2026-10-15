"""Write design.csv + a command file for virtual-ecologist runs with REEIM_vesurvey.slim (a copy of REEIM.slim that also
writes survey_*.txt: 2000 random females every 5 years 1900-2020 with age, length, maturity and first-time maturing).
Usage: python3 make_vesurvey.py <out dir> <base dir holding REEIM/input and REEIM_vesurvey.slim> <cmds file name> <first rep> 'H,SLT,w,gen,label' ...
gen: default | more (Ve 0.02, sd 0.1) | extreme (Ve 0.01, sd 0.1, u 1e-5), as in make_stress.py. Appends to design.csv.
Analyse with analyse_vesurvey.R."""
import sys, os
V, B, cname, first = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
gen = {"default": (0.05, 0.05, 1e-6), "more": (0.02, 0.1, 1e-6), "extreme": (0.01, 0.1, 1e-5)}
os.makedirs(os.path.join(V, "out/plaice"), exist_ok=True)
dp = os.path.join(V, "design.csv"); new = not os.path.exists(dp)
with open(dp, "a") as d, open(os.path.join(V, cname), "w") as c:
    if new: d.write("set,rep,H,f,w,Shmax,SLT,gen,Ve,sd,u,label\n")
    for i, spec in enumerate(sys.argv[5:], first):
        H, SLT, w, g, lab = spec.split(",", 4); Ve, sd, u = gen[g]
        d.write(f"plaice,{i},{H},0.1,{w},117,{SLT},{g},{Ve},{sd},{u},{lab}\n")
        c.write(f"slim -s {5000+i} -d \"folder='{V}/out/plaice'\" -d \"inputDir='{B}/REEIM/input/plaice'\" -d Tdelta=10000 -d rep={i} "
                f"-d H={H} -d f=0.1 -d w={w} -d Shmax=117 -d SLT={SLT} -d Ve={Ve} -d sd={sd} -d u={u} {B}/REEIM_vesurvey.slim > /dev/null 2>&1\n")
