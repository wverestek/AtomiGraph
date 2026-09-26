# Generate two LAMMPS data files (atom_style full) for a united-atom polyethylene chain:
# frame 0: intact chain; frame 1000: two well-separated chain scissions.
import math, sys
N = 30                          # united atoms: CH3-(CH2)28-CH3
BROKEN = {(8, 9), (21, 22)}     # scission sites, > 2*rxn_bond_cutoff bonds apart
b, theta = 1.54, math.radians(112.0)
dx, dy = b*math.sin(theta/2), b*math.cos(theta/2)
out = sys.argv[1]

def write(fname, ts, bonds):
    with open(fname, "w") as f:
        f.write(f"LAMMPS data file, united-atom PE chain, timestep = {ts}\n\n")
        f.write(f"{N} atoms\n2 atom types\n{len(bonds)} bonds\n1 bond types\n\n")
        f.write("0.0 50.0 xlo xhi\n-10.0 10.0 ylo yhi\n-10.0 10.0 zlo zhi\n\n")
        f.write("Masses\n\n1 14.027\n2 15.035\n\n")
        f.write("Atoms # full\n\n")
        for i in range(1, N+1):
            t = 2 if i in (1, N) else 1           # 1 = CH2, 2 = CH3 (chain ends)
            x, y = 2.0 + (i-1)*dx, (i % 2)*dy
            f.write(f"{i} 1 {t} 0.0 {x:.4f} {y:.4f} 0.0000 0 0 0\n")
        f.write("\nBonds\n\n")
        for k, (i, j) in enumerate(bonds, 1):
            f.write(f"{k} 1 {i} {j}\n")

chain = [(i, i+1) for i in range(1, N)]
write(f"{out}/pe_chain.0.data", 0, chain)
write(f"{out}/pe_chain.1000.data", 1000, [p for p in chain if p not in BROKEN])
