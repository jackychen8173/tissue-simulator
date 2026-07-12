import math

# Full PIN + PD auxin canalization model (Holloway, Eiriksson & Wenzel 2025,
# Front. Plant Sci. 16:1621815), implemented on the semicircular leaf tissue.
#
# Cell variables:
#   0 = volume tracker (for Division::ShortestPath2D)
#   1 = type (0 = normal, 1 = source/outer margin ring, 2 = sink/base ring)
#   2 = Aprec (auxin precursor; ramps in source cells, drives auxin production)
#   3 = A (auxin)
#   4 = P (PIN, cytoplasmic)
#
# Wall variables (index 0 is always the built-in rest length):
#   1, 2 = Pij, Pji (PIN membrane pair, directional)
#   3    = Dij (plasmodesmata cross-sectional area, symmetric; only this index
#          is touched by any reaction)
#   4    = unused mirror slot, present only so the VTU writer's paired-output
#          convention (it always reads wall variables two at a time) doesn't
#          read one index out of bounds with an odd variable count
#   5, 6 = saved PIN-transport flux pair (directional), feeds the WTF term

def createSemicircle(nRings, nCols, outputName):

    R_max = 10.0
    radii = [R_max * (r+1) / nRings for r in range(nRings)]
    angles = [i * math.pi / nCols for i in range(nCols + 1)]

    verts = []
    # 1. Create the base vertices evenly spaced along the FLAT top axis (y = 0)
    # This prevents MoveVertexRadially from warping a curved inner arc into the flat zone.
    for i in range(nCols + 1):
        x_flat = - (radii[0] * 0.1) + i * ((2 * radii[0] * 0.1) / nCols)
        verts.append((x_flat, 0.0))

    # 2. Create the concentric rings (True arcs)
    for r in range(nRings):
        R = radii[r]
        for i in range(nCols + 1):
            a = angles[i]
            verts.append((-R * math.cos(a), R * math.sin(a)))

    nVerts = len(verts)
    nCells = nRings * nCols

    def base_vert(i):   return i
    def arc_vert(r, i): return (r+1)*(nCols+1) + i
    def cell_idx(r, c): return r * nCols + c

    walls = []
    for c in range(nCols):
        walls.append((-1, cell_idx(0,c), base_vert(c), base_vert(c+1)))
    for r in range(nRings-1):
        for c in range(nCols):
            walls.append((cell_idx(r,c), cell_idx(r+1,c), arc_vert(r,c), arc_vert(r,c+1)))
    for c in range(nCols):
        r = nRings-1
        walls.append((-1, cell_idx(r,c), arc_vert(r,c), arc_vert(r,c+1)))
    for r in range(nRings):
        v_inner = base_vert(0) if r==0 else arc_vert(r-1,0)
        walls.append((-1, cell_idx(r,0), v_inner, arc_vert(r,0)))
    for r in range(nRings):
        for c in range(1, nCols):
            v_inner = base_vert(c) if r==0 else arc_vert(r-1,c)
            walls.append((cell_idx(r,c-1), cell_idx(r,c), v_inner, arc_vert(r,c)))
    for r in range(nRings):
        v_inner = base_vert(nCols) if r==0 else arc_vert(r-1,nCols)
        walls.append((-1, cell_idx(r,nCols-1), v_inner, arc_vert(r,nCols)))

    nWalls = len(walls)

    # Initial wall values: PIN pair = 0, Dij seeded nonzero (PD transport/feedback
    # is inert while Dij<=0, see DiffusionConductiveSimple), flux-save pair = 0.
    Dij_init = 0.5

    # Initial cell values
    volume_init = 0.9
    Aprec_source_init = 0.5   # matches paper's initial Z1 precursor spark
    PIN_cytoplasmic_init = 0.1

    with open(outputName, 'w') as f:
        f.write(f"{nCells} {nWalls} {nVerts}\n")
        for idx, (c1,c2,v1,v2) in enumerate(walls):
            f.write(f" {idx} {c1} {c2} {v1} {v2}\n")
        f.write(" \n")
        f.write(f" {nVerts} 2\n")
        for (x,y) in verts:
            f.write(f"{x:.6f} {y:.6f}\n")
        # 6 extra wall columns beyond the built-in length: Pij Pji Dij (mirror) flux1 flux2
        f.write(f"\n{nWalls} 1 6\n ")
        for (c1,c2,v1,v2) in walls:
            vx1,vy1 = verts[v1]
            vx2,vy2 = verts[v2]
            rest = math.sqrt((vx2-vx1)**2+(vy2-vy1)**2)
            f.write(f"{rest:.6f} 0.0 0.0 {Dij_init} {Dij_init} 0.0 0.0\n")
        f.write(f"\n{nCells} 5\n")
        for r in range(nRings):
            cellType = 1 if r == nRings - 1 else (2 if r == 0 else 0)
            for c in range(nCols):
                aprec = Aprec_source_init if cellType == 1 else 0.0
                f.write(f"{volume_init} {cellType} {aprec} 0.0 {PIN_cytoplasmic_init}\n")

    print(f"Written {outputName}: {nCells} cells, {nWalls} walls, {nVerts} vertices")


createSemicircle(
    nRings=6, nCols=20,
    outputName="semicircle_leaf2.init"
)


# Reaction parameters — stepped toward Holloway et al. (2025) paper values.
# Key changes from previous run (4x14, t=100):
#   - Tissue: 6 rings x 20 cols (120 cells, up from 56) for longer canal. path
#   - T_pin: 1.0 -> 3.0  (paper: 6; stronger PIN-driven transport)
#   - D_pd:  0.05 -> 0.15 (paper: 0.8; more PD diffusion while T still > D)
#   - k_U:   0.05 -> 0.1  (paper: 4e-3 rescaled; stronger up-gradient)
#   - k_Wl:  0.1  -> 0.2  (paper: 3e-3 rescaled; stronger with-the-flux)
#   - k_Wq:  0.5  -> 0.15 (reduce — was destabilizing at high T)
#   - k_off: 0.02 -> 0.05 (more PIN redistribution, less early locking)
#   - pinpr: 0.02 -> 0.05 (more cytoplasmic PIN available where auxin is high)
#   - gamma_over_alpha: 2.5 -> 5.0 (tighter Dij ceiling, narrower tracks)
#   - beta_pd: 0.01 -> 0.02 (more background PD to keep off-track walls open)
#   - V_threshold: 8.0 -> 6.0 (smaller initial cell area with 6 rings)
#   - Run time: t=100 -> t=200
params = {
    "k_growth_radial": 0.006,
    "k_growth_wall":   0.048,
    "K_spring":        0.01,
    "K_adh":           1.0,
    # Eq1: auxin
    "auxpr":       0.05,   # production rate from Aprec
    "auxdec":      0.01,   # background decay
    "auxdec_sink": 0.15,   # extra decay at sink ring
    "aprec_ramp":  0.02,   # precursor ramp rate in source cells
    "T_pin":       3.0,    # PIN transport permeability (paper: 6)
    # Eq2: plasmodesmata (via DiffusionConductiveSimple, p3=0)
    "D_pd":            0.15,  # PD passive diffusion rate (paper: 0.8)
    "alpha_pd":        0.02,  # flux-feedback growth rate
    "gamma_over_alpha": 5.0,  # gamma/alpha ratio -> sets Dij ceiling
    "beta_pd":         0.02,  # background PD production on all walls
    # Eq3: PIN cytoplasmic
    "pinpr":  0.05,
    "pindec": 0.02,
    # Eq4: PIN allocation (UTG + WTF)
    "k_U":   0.1,
    "k_off": 0.05,
    "K_utg": 2.0,
    "k_Wq":  0.15,
    "k_Wl":  0.2,
    # Division / removal
    "V_threshold": 6.0,
    "R_removal":   500.0,
}

model = f"""\
14
2
0
MoveVertexRadially
2 0
{params['k_growth_radial']}
1
WallGrowth::Stress
4 2 1 1
{params['k_growth_wall']}
0.0
1
1
0
1
VertexFromWallSpring
2 1 1
{params['K_spring']}
{params['K_adh']}
0
CenterCOM 0 0
Creation::FromType
2 2 1 1
{params['aprec_ramp']}
1
2
1
Creation::One
1 2 1 1
{params['auxpr']}
3
2
Degradation::One
1 1 1
{params['auxdec']}
3
Degradation::FromType
2 2 1 1
{params['auxdec_sink']}
2
3
1
PINSaturatingTransport
1 3 1 1 1
{params['T_pin']}
3
1
5
Diffusion::ConductiveSimple
5 2 1 2
{params['D_pd']}
{params['alpha_pd']}
2
0
{params['gamma_over_alpha']}
3
3
4
Creation::OneWall
1 1 2
{params['beta_pd']}
3
4
Creation::One
1 2 1 1
{params['pinpr']}
4
3
Degradation::One
1 1 1
{params['pindec']}
4
MembraneCycling::UTGWTF
5 3 2 1 1
{params['k_U']}
{params['k_off']}
{params['K_utg']}
{params['k_Wq']}
{params['k_Wl']}
3
4
1
5
Division::ShortestPath2D
4 1 1
{params['V_threshold']}
1.0
0.1
1
0
RemovalOutsideRadius
1 0
{params['R_removal']}
"""

with open("semicircle_leaf2.model", "w") as f:
    f.write(model)
print("Written semicircle_leaf2.model")


solver = """\
RK5Adaptive
0 112
2 56
0.5 1e-5
"""

with open("semicircle_leaf2.solver", "w") as f:
    f.write(solver)
print("Written semicircle_leaf2.solver")
