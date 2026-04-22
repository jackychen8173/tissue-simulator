# Leaf Vein Simulation — Semicircular Monocot Leaf Model

A computational simulation of auxin transport and PIN-mediated polar auxin transport
in a semicircular monocot leaf geometry, implemented using the
[Tissue simulator](https://gitlab.com/slcu/teamhj/tissue) (Jonsson lab, SLCU Cambridge).

This project implements and extends the model described in:

> Holloway DM, Eiriksson TK, Wenzel CL (2025). "The role of auxin transport through
> plasmodesmata in leaf vein canalization and patterning."
> *Frontiers in Plant Science* 16:1621815.

---

## Overview

The simulation models how the plant hormone auxin creates narrow provascular strands
(future leaf veins) in a growing semicircular leaf tissue. Auxin is produced at the
apex and transported toward the base via PIN protein-mediated polar auxin transport,
canalizing into a narrow stripe that marks the future vein location.

---

## Repository Structure
```
├── src/                          # Tissue simulator source code (MIT License, SLCU/teamHJ)
├── bin/                          # Compiled simulator binary (after building)
├── examples/
│   └── tutorials/
│       └── semicircle_leaf/      # Our custom leaf model
│           ├── generate_tissue_files.py   # Generates all 3 input files
│           ├── semicircle_leaf.model      # Reaction rules
│           ├── semicircle_leaf.init       # Starting geometry
│           ├── semicircle_leaf.solver     # Time and output settings
│           └── vtk/                       # Simulation output (generated on run)
├── LICENSE                       # MIT License (original Tissue simulator)
└── README.md
```


---

## Requirements

- Linux or WSL (Windows Subsystem for Linux)
- C++ build tools
- Boost C++ libraries
- Python 3
- [Paraview](https://www.paraview.org/) for visualization

---

## Installation

### 1. Install WSL (Windows only)

Open PowerShell as Administrator:
```bash
wsl --install
```
Restart when prompted. Then open Ubuntu from the Start menu.

### 2. Install Dependencies

Inside WSL:
```bash
sudo apt update
sudo apt install build-essential cmake git python3 python3-pip
sudo apt install libboost-all-dev
```

### 3. Clone This Repository

```bash
git clone https://github.com/YOUR_USERNAME/leaf-vein-simulation.git
cd leaf-vein-simulation
```

### 4. Build the Tissue Simulator

```bash
cd src
make
cd ..
```

The compiled binary will be at `bin/simulator`.

### 5. Install Paraview

Download and install Paraview from [paraview.org](https://www.paraview.org/).

---

## Running the Simulation

### Step 1 — Generate Input Files

```bash
cd examples/tutorials/semicircle_leaf
python3 generate_tissue_files.py
```

This creates three files:
- `semicircle_leaf.init` — semicircular leaf geometry (56 cells, 4 rings, 14 columns)
- `semicircle_leaf.model` — auxin/PIN reaction rules using `AuxinModel1S`
- `semicircle_leaf.solver` — RK5 adaptive solver settings

### Step 2 — Run the Simulation

```bash
../../bin/simulator semicircle_leaf.model semicircle_leaf.init semicircle_leaf.solver
```

Wait for `Simulation done.` Output files are written to `vtk/`.

### Step 3 — Visualize in Paraview

1. Open Paraview
2. File → Open → navigate to `vtk/` → open `tissue.pvd`
3. Click **Apply**
4. In the color dropdown, select:
   - `cell variable 4` — auxin concentration
   - `cell variable 5` — PIN cytoplasm
   - `wall variable 1` — PIN membrane (polarization)
5. Press Play to animate through time steps

---

## Model Description

### Geometry

The semicircular leaf has a flat base (y=0) and an arc curving downward.
Cells are arranged in 4 radial rings × 14 angular columns = 56 cells total.
Auxin is seeded at high concentration in the 2 central apex cells.

### Cell Variables

| Paraview Name | Index | Description |
|---|---|---|
| cell variable 4 | 0 (user) | Auxin concentration Ai |
| cell variable 5 | 1 (user) | PIN cytoplasm pool Pi |
| cell variable 6 | 2 (user) | X boundary marker (unused) |
| cell variable 7 | 3 (user) | M membrane marker |

### Wall Variables

| Paraview Name | Description |
|---|---|
| wall length | Mechanical rest length |
| wall variable 1 | PIN membrane side 1 (Pij) |
| wall variable 2 | PIN membrane side 2 (Pji) |

### Reaction Rules

**`AuxinModel1S`** — combined auxin transport and PIN cycling rule implementing:
- Auxin production (`c_A = 0.005`) and decay (`d_A = 0.005`)
- PIN polarization: 90% polar (`k_p = 0.9`), feedback strength `f_p = 0.3`
- Active PIN transport rate `T = 1.3`, passive diffusion `D = 0.002`
- PIN production (`c_P = 0.005`) and decay (`d_P = 0.005`)

**`WallMechanics::Spring`** — wall elasticity maintaining tissue shape

**`CenterCOM`** — keeps tissue centred at origin

---

## Customization

Edit parameters in `generate_tissue_files.py`:

| Parameter | Location | Effect |
|---|---|---|
| `nRings` | `createSemicircle()` call | More rings = finer radial resolution |
| `nCols` | `createSemicircle()` call | More columns = finer angular resolution |
| `R_max` | inside `createSemicircle()` | Starting tissue size |
| `T` | model string | Active transport rate — increase for stronger canalization |
| `D` | model string | Passive diffusion — decrease for narrower auxin stripe |
| `k_p` | model string | Fraction of polarized PIN — 0.9 = 90% polar |
| End time | solver string | Longer simulation = more developed pattern |

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `boost/mpl/fold.hpp not found` | `sudo apt install libboost-all-dev` |
| `Cannot open file` | Run simulator from inside the `semicircle_leaf/` folder |
| No `.vtu` files produced | Check solver output format is `2` |
| Auxin not spreading | Check variable indices: auxin=4, PIN=5 in Paraview |
| Simulation unstable | Reduce `T` and `c_A` in model string |

---

## License

The Tissue simulator source code is licensed under the
[MIT License](LICENSE) © 2017 SLCU / teamHJ.

The semicircular leaf model and `generate_tissue_files.py` were developed as part of
a graduate research project at BCIT implementing the Holloway et al. (2025) auxin
transport model in a monocot leaf geometry.

---

## References

- Holloway DM, Eiriksson TK, Wenzel CL (2025). The role of auxin transport through
  plasmodesmata in leaf vein canalization and patterning.
  *Frontiers in Plant Science* 16:1621815.
- Caggiano MP, Yu X, Bhatia N, Larsson A, et al. (2017). Cell type boundaries organize
  plant development. *eLife* 6:e27421.
- Merks RMH, Guravage M, Inzé D, Beemster GTS (2011). VirtualLeaf: An open-source
  framework for cell-based modeling of plant tissue growth and development.
  *Plant Physiology* 155:656–666.