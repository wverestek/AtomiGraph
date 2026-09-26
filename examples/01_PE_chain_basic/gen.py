# Generate LAMMPS data files (atom_style full) for a united-atom polyethylene chain.
# Each frame changes the bond topology of the previous frame by one scenario:
#   timestep    0: intact chain 1-2-...-40
#   timestep 1000: scission     8-9 and 21-22 (well separated -> 2 reactions)
#   timestep 2000: formation    8-22          (recombination of two chain ends -> 1 reaction)
#   timestep 3000: flip         30-31 -> 31-15 (bond broken and new bond formed -> 1 reaction)
#   timestep 4000: flip 33-34 -> 34-12 plus scission 36-37, 2 bonds apart
#                  (cores separate: 2 reactions for cutoff 0, merged: 1 reaction for cutoff >= 1)
# Graph test data: branch points keep the CH2 type, geometry is not updated.
import math, sys
N = 40                          # united atoms: CH3-(CH2)38-CH3
b, theta = 1.54, math.radians(112.0)
dx, dy = b*math.sin(theta/2), b*math.cos(theta/2)
out = sys.argv[1] if len(sys.argv) > 1 else "."

# (timestep, bonds broken, bonds formed) relative to the previous frame
STEPS = [
    (1000, [(8, 9), (21, 22)], []),
    (2000, [], [(8, 22)]),
    (3000, [(30, 31)], [(31, 15)]),
    (4000, [(33, 34), (36, 37)], [(34, 12)]),
]

def write(fname, ts, bonds):
    with open(fname, "w") as f:
        f.write(f"LAMMPS data file, united-atom PE chain, timestep = {ts}\n\n")
        f.write(f"{N} atoms\n2 atom types\n{len(bonds)} bonds\n1 bond types\n\n")
        f.write("0.0 70.0 xlo xhi\n-10.0 10.0 ylo yhi\n-10.0 10.0 zlo zhi\n\n")
        f.write("Masses\n\n1 14.027\n2 15.035\n\n")
        f.write("Atoms # full\n\n")
        for i in range(1, N+1):
            t = 2 if i in (1, N) else 1           # 1 = CH2, 2 = CH3 (chain ends)
            x, y = 2.0 + (i-1)*dx, (i % 2)*dy
            f.write(f"{i} 1 {t} 0.0 {x:.4f} {y:.4f} 0.0000 0 0 0\n")
        f.write("\nBonds\n\n")
        for k, (i, j) in enumerate(bonds, 1):
            f.write(f"{k} 1 {i} {j}\n")

bonds = [(i, i+1) for i in range(1, N)]
write(f"{out}/pe_chain.0.data", 0, bonds)
for ts, broken, formed in STEPS:
    bonds = [p for p in bonds if p not in broken] + formed
    write(f"{out}/pe_chain.{ts}.data", ts, bonds)
