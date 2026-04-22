import math
import random

def createSquareGrid(nGrid, outputName, numCellVariables, cellsMean,
    cellsSpread, numWallVariables, wallsMean, wallsSpread):
    n = nGrid
    nCells = n**2
    nWALLS = 4 + 6*(n-1) + 2*(n-1)**2
    nVERTEX = (n+1)**2
    c = 0
    ncount = n-1
    vertexone = 0
    vertextwo = 1
    chan = open(outputName, "w")
    chan.write(str(nCells)+" "+str(nWALLS)+" "+str(nVERTEX)+"\n")

    while c <= ncount:
        chan.write(" "+str(c)+" "+str(c)+" -1 "+str(vertexone)+" "+str(vertextwo)+"\n")
        c += 1
        vertexone += 1
        vertextwo += 1

    cellone = -1
    celltwo = 0
    vertexone = n+1
    vertextwo = 0
    wallcounter = 1
    celloneA = 0
    celltwoA = 0
    vertexoneA = n+1
    vertextwoA = n+2

    while wallcounter <= n-1:
        ncount = ncount+n+1
        chan.write(" "+str(c)+" "+str(celltwo)+" -1 "+str(vertextwo)+" "+str(vertexone)+"\n")
        c += 1
        vertexone += 1
        vertextwo += 1
        celltwo += 1
        cellone += 1
        while c <= ncount-1:
            chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
            c += 1
            vertexone += 1
            vertextwo += 1
            celltwo += 1
            cellone += 1
        chan.write(" "+str(c)+" "+str(cellone)+" -1 "+str(vertextwo)+" "+str(vertexone)+"\n")
        vertexone += 1
        vertextwo += 1
        ncount += n
        c += 1
        celltwoA = celltwo
        while c <= ncount:
            chan.write(" "+str(c)+" "+str(celloneA)+" "+str(celltwoA)+" "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
            c += 1
            vertextwoA += 1
            vertexoneA += 1
            celltwoA += 1
            celloneA += 1
        vertextwoA += 1
        vertexoneA += 1
        wallcounter += 1

    ncount = ncount+n+1
    chan.write(" "+str(c)+" "+str(celltwo)+" -1 "+str(vertextwo)+" "+str(vertexone)+"\n")
    c += 1
    vertexone += 1
    vertextwo += 1
    celltwo += 1
    cellone += 1
    while c <= ncount-1:
        chan.write(" "+str(c)+" "+str(cellone)+" "+str(celltwo)+" "+str(vertextwo)+" "+str(vertexone)+"\n")
        c += 1
        vertexone += 1
        vertextwo += 1
        celltwo += 1
        cellone += 1
    chan.write(" "+str(c)+" "+str(cellone)+" -1 "+str(vertextwo)+" "+str(vertexone)+"\n")
    vertexone += 1
    vertextwo += 1
    ncount += n
    c += 1
    while c <= nWALLS-1:
        chan.write(" "+str(c)+" "+str(celloneA)+" -1 "+str(vertexoneA)+" "+str(vertextwoA)+"\n")
        c += 1
        vertextwoA += 1
        vertexoneA += 1
        celltwoA += 1
        celloneA += 1

    chan.write(" \n")
    chan.write(" "+str(nVERTEX)+" 2\n")
    ycoord = 0
    xcoord = 0
    while xcoord <= n:
        while ycoord <= n:
            chan.write(str(xcoord)+" "+str(ycoord)+"\n")
            ycoord += 1
        ycoord = 0
        xcoord += 1

    wallLength = 1.0
    chan.write("\n"+str(nWALLS)+" 1 "+str(numWallVariables*2)+"\n ")
    for iw in range(nWALLS):
        chan.write(str(wallLength)+" ")
        for j in range(numWallVariables):
            d1 = wallsMean[j]+wallsSpread[j]*(2*random.uniform(0,1)-1)
            d2 = wallsMean[j]+wallsSpread[j]*(2*random.uniform(0,1)-1)
            chan.write(str(d1)+" "+str(d2)+" ")
        chan.write("\n")

    chan.write("\n"+str(nCells)+" "+str(numCellVariables+4)+"\n")
    # seed auxin in centre cells
    centre = n // 2
    for r in range(n):
        for col in range(n):
            auxin = 0.5 if (r == centre and col in (centre-1, centre)) else 0.0
            chan.write("1 0 1 0.9 "+str(auxin)+" ")
            for j in range(numCellVariables-1):
                d = cellsMean[j]+cellsSpread[j]*(2*random.uniform(0,1)-1)
                chan.write(str(d)+" ")
            chan.write("\n")
    chan.close()
    print(f"Written {outputName}: {nCells} cells, {nWALLS} walls, {nVERTEX} vertices")


createSquareGrid(
    nGrid=8,
    outputName="square_leaf.init",
    numCellVariables=5,
    cellsMean=[0.0, 0.1, 0.0, 0.0, 0.0],
    cellsSpread=[0.0]*5,
    numWallVariables=2,
    wallsMean=[0.0, 0.0],
    wallsSpread=[0.0, 0.0]
)


model = """\
10
0
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
3

WallMechanics::Spring
2 1 1
0.01
1.0
0

CenterCOM 0 0

Creation::Zero
1 1 1
0.02
0

Degradation::One
1 1 1
0.05
0

DiffusionActiveTransportCell
2 2 1 1
0.1
2.0
0
1

Creation::One
1 2 1 1
0.05
1
0

Degradation::One
1 1 1
0.02
1

MembraneCycling::CellUpTheGradientNonLinear
4 2 2 1
0.1
0.05
1.0
2.0
0
1
1
"""

with open("square_leaf.model", "w") as f:
    f.write(model)
print("Written square_leaf.model")


solver = """\
RK5Adaptive
0 3000
2
100
0.5 1e-3
"""

with open("square_leaf.solver", "w") as f:
    f.write(solver)
print("Written square_leaf.solver")