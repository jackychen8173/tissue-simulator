import math
import random

def createSemicircle(nRings, nCols, outputName, numCellVariables,
    cellsMean, cellsSpread, numWallVariables, wallsMean, wallsSpread):

    R_max = 10.0
    radii = [R_max * (r+1) / nRings for r in range(nRings)]
    angles = [i * math.pi / nCols for i in range(nCols + 1)]

    R_inner = R_max * 0.05  # small inner arc to avoid zero-length walls
    verts = []
    for i in range(nCols + 1):
        a = angles[i]
        verts.append((-R_inner * math.cos(a), -R_inner * math.sin(a)))

    # 2. Create the concentric rings
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
        for r in range(nRings):
            for c in range(nCols):
                extra = ""
                for j in range(numCellVariables):
                    val = cellsMean[j] + cellsSpread[j] * (2*random.uniform(0,1)-1)
                    extra += f" {val}"
                f.write(f"1 0 1 0.9{extra}\n")

    print(f"Written {outputName}: {nCells} cells, {nWalls} walls, {nVerts} vertices")


createSemicircle(
    nRings=4, nCols=14,
    outputName="semicircle_leaf2.init",
    numCellVariables=5,
    cellsMean=[0.0]*5,
    cellsSpread=[0.0]*5,
    numWallVariables=2,
    wallsMean=[0.0, 0.0],
    wallsSpread=[0.0, 0.0]
)


model = """\
3
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
Division::ShortestPath2D
4 1 1
5.0
1.0
0.1
1
3
RemovalOutsideRadius
1 0
15.0
"""

with open("semicircle_leaf2.model", "w") as f:
    f.write(model)
print("Written semicircle_leaf2.model")


solver = """\
RK5Adaptive
0 5000
2 100
0.5 1e-5
"""

with open("semicircle_leaf2.solver", "w") as f:
    f.write(solver)
print("Written semicircle_leaf2.solver")