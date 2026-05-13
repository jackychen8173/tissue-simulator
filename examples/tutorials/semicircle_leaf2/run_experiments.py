import subprocess
import os
import shutil
import math
import random

def createSemicircle(nRings, nCols, outputName, numCellVariables,
    cellsMean, cellsSpread, numWallVariables, wallsMean, wallsSpread, auxin_apex=50.0):

    R_max = 200.0
    radii = [R_max * (r+1) / nRings for r in range(nRings)]
    angles = [i * math.pi / nCols for i in range(nCols + 1)]

    verts = []
    # This replaces the flat x = -R_max to R_max loop
    for i in range(nCols + 1):
        a = angles[i]
        # Calculate the starting inner radius (r=0 would be 50.0). 
        # Using a very small radius for the origin avoids distortion.
        r_inner = radii[0] * 0.1  
        verts.append((-r_inner * math.cos(a), -r_inner * math.sin(a)))

    # 2. Create the rest of the rings 
    for r in range(nRings):
        R = radii[r]
        for i in range(nCols + 1):
            a = angles[i]
            verts.append((-R * math.cos(a), -R * math.sin(a)))

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

    with open(outputName, 'w') as f:
        f.write(f"{nCells} {nWalls} {nVerts}\n")
        for idx, (c1,c2,v1,v2) in enumerate(walls):
            f.write(f" {idx} {c1} {c2} {v1} {v2}\n")
        f.write(" \n")
        f.write(f" {nVerts} 2\n")
        for (x,y) in verts:
            f.write(f"{x:.6f} {y:.6f}\n")
        f.write(f"\n{nWalls} 1 {2*numWallVariables}\n ")
        for (c1,c2,v1,v2) in walls:
            vx1,vy1 = verts[v1]
            vx2,vy2 = verts[v2]
            rest = math.sqrt((vx2-vx1)**2+(vy2-vy1)**2)
            f.write(f"{rest:.6f} ")
            for j in range(numWallVariables):
                d1 = wallsMean[j]+wallsSpread[j]*(2*random.uniform(0,1)-1)
                d2 = wallsMean[j]+wallsSpread[j]*(2*random.uniform(0,1)-1)
                f.write(f"{d1} {d2} ")
            f.write("\n")
        f.write(f"\n{nCells} {numCellVariables + 4}\n")
        apex_cols = set(range(nCols//2 - 1, nCols//2 + 1))
        for r in range(nRings):
            for c in range(nCols):
                auxin = auxin_apex if (r == nRings-1 and c in apex_cols) else 0.0
                f.write(f"1 0 1 0.9 {auxin} 0.1 0.0 0.0 0.0\n")

    print(f"Written {outputName}: {nCells} cells, {nWalls} walls, {nVerts} vertices")


# Parameter combinations to test
experiments = [
    {"name": "low_transport",  "c_A": 0.05, "d_A": 0.001, "T": 0.5, "D": 0.002, "kU": 0.3, "nRings": 4, "nCols": 14},
    {"name": "med_transport",  "c_A": 0.05, "d_A": 0.001, "T": 1.3, "D": 0.002, "kU": 0.3, "nRings": 4, "nCols": 14},
    {"name": "high_transport", "c_A": 0.05, "d_A": 0.001, "T": 6.0, "D": 0.002, "kU": 0.3, "nRings": 6, "nCols": 20},
    {"name": "paper_values",   "c_A": 0.03, "d_A": 0.05,  "T": 6.0, "D": 0.002, "kU": 4e-3, "nRings": 8, "nCols": 28},
]

SIMULATOR = "/mnt/c/Users/jchen/tissue-simulator/bin/simulator"
BASE_DIR = "/mnt/c/Users/jchen/tissue-simulator/examples/tutorials/semicircle_leaf2"

for exp in experiments:
    print(f"\nRunning experiment: {exp['name']}")

    # Generate init file
    createSemicircle(
        nRings=exp['nRings'], nCols=exp['nCols'],
        outputName=f"{BASE_DIR}/semicircle_leaf2.init",
        numCellVariables=5,
        cellsMean=[0.0]*5,
        cellsSpread=[0.0]*5,
        numWallVariables=2,
        wallsMean=[0.0, 0.0],
        wallsSpread=[0.0, 0.0],
        auxin_apex=50.0
    )

    # Generate model file
    model = f"""\
9
2
0
MoveVertexRadially
2 0
0.001
1
WallGrowth::Stress
4 2 1 1
0.005
0.0
1
1
0
1
VertexFromWallSpring
2 1 1
0.01
1.0
0
CenterCOM 0 0
Creation::Zero
1 1 1
{exp['c_A']}
4
Degradation::One
1 1 1
{exp['d_A']}
4
DiffusionActiveTransportCell
2 2 1 1
{exp['D']}
{exp['T']}
4
1
Creation::One
1 2 1 1
0.005
5
4
Degradation::One
1 1 1
0.005
5
MembraneCycling::CellUpTheGradientNonLinear
4 2 2 1
{exp['kU']}
0.1
1.0
2.0
4
5
1
Division::ShortestPath2D
4 1 1
1.0
1.0
0.3
1
3
RemovalOutsideRadius
1 0
220.0
"""

    with open(f"{BASE_DIR}/semicircle_leaf2.model", "w") as f:
        f.write(model)

    # Generate solver file
    solver = """\
RK5Adaptive
0 10000
2 50
0.05 1e-5
"""
    with open(f"{BASE_DIR}/semicircle_leaf2.solver", "w") as f:
        f.write(solver)

    # Run simulator
    result = subprocess.run(
        [SIMULATOR, "semicircle_leaf2.model", "semicircle_leaf2.init", "semicircle_leaf2.solver"],
        cwd=BASE_DIR,
        capture_output=True,
        text=True
    )
    print(result.stdout[-200:] if result.stdout else "No output")

    # Save vtk output to named folder
    output_dir = f"{BASE_DIR}/results/{exp['name']}"
    os.makedirs(output_dir, exist_ok=True)
    vtk_src = f"{BASE_DIR}/vtk"
    if os.path.exists(vtk_src):
        shutil.copytree(vtk_src, f"{output_dir}/vtk", dirs_exist_ok=True)
        print(f"Results saved to results/{exp['name']}/vtk/")

print("\nAll experiments complete.")