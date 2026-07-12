# Project Changes Log

Summary of every modification made to the original VirtualLeaf codebase for this project (implementing Holloway, Eiriksson & Wenzel 2025, *Front. Plant Sci.* 16:1621815).

---

## ParaView Variable Guide

Open `results/<run>/tissue.pvd` in ParaView to load the full time series.

### Cells (`VTK_cells*.vtu`)

The cell data writer in `src/VTUostream.cc` always packs cell variables in this fixed layout:

| ParaView name | Internal variables | What it shows |
|---|---|---|
| **cell vector** | cell[0], cell[1], cell[2] | 3-component vector: [volume tracker, cell type, Aprec (auxin precursor)] — mainly useful for checking cell types (component Y = 0 normal, 1 source/outer ring, 2 sink/base ring) |
| **cell vector length** | cell[3] | **Auxin (A)** — the main signaling molecule; use this to see the auxin gradient and canal |
| **cell variable 4** | cell[4] | **Cytoplasmic PIN (P)** — mobile pool of PIN protein inside the cell; higher where auxin is high |

> **To visualize the auxin pattern:** color by **cell vector length**
> **To visualize cytoplasmic PIN:** color by **cell variable 4**
> **To identify cell types:** color by **cell vector**, component index 1 (Y component)

---

### Walls (`VTK_walls*.vtu`)

The wall writer uses a **paired** format: consecutive wall variable slots are grouped into a 2-component array. So each "wall variable N" in ParaView holds two values — one for each direction across the wall (i→j and j→i).

| ParaView name | Internal variables | What it shows |
|---|---|---|
| **wall length** | built-in rest length | Current resting length of each wall; grows over time via `WallGrowth::Stress` |
| **wall variable 0** | [wall[0], wall[1]] = [Pij, Pji] | **PIN on each face of the wall** — component 0 = PIN pointing from cell1 to cell2, component 1 = PIN pointing from cell2 to cell1. Asymmetry between components shows PIN polarity. |
| **wall variable 1** | [wall[2], wall[3]] = [Dij, Dij_mirror] | **Plasmodesmata area (Dij)** — both components are now equal (Dij = Dji, symmetric pores). Larger values = wider PD channels = more diffusion. Walls in canalized tracks accumulate Dij over time. |
| **wall variable 2** | [wall[4], wall[5]] = [flux_ij, flux_ji] | **Saved PIN-transport flux** — magnitude of auxin pumped across the wall each step, one component per direction. High values = walls carrying active PIN-mediated auxin flow. |

> **To see PIN polarity:** color walls by **wall variable 0**, then split components to compare Pij vs Pji on each wall
> **To see PD canalization:** color walls by **wall variable 1** (both components now equal) — brighter walls are wider PD channels
> **To see which walls carry active flux:** color walls by **wall variable 2**

---

### Why the wall variables are paired

`VTUostream.cc` (line 1306) iterates wall slots two at a time (`i += 2`) and names each output array `"wall variable " << i/2`. This halves the number of visible arrays in ParaView. The pairing was designed for directional quantities that naturally have an i→j and j→i value — PIN allocation and PD area both fit this structure.

---

## 1. New C++ Reaction Classes

Five new reaction classes were written and added to the existing source files. They follow the exact same conventions as the existing VirtualLeaf classes (same `derivs()` signature, same `variableIndex()` / `parameter()` API, same conservation pattern for wall↔cell exchange).

### `Creation::FromType`
**File:** [src/creation.cc](src/creation.cc), [src/creation.h](src/creation.h)

Produces a molecule only in cells whose "type" variable matches a specified value. Used to ramp up auxin precursor (`Aprec`) exclusively in the outer source ring.

Modeled after `Creation::One`; added a type-check condition inside the cell loop.

---

### `Creation::OneWall`
**File:** [src/creation.cc](src/creation.cc), [src/creation.h](src/creation.h)

Produces a molecule at a constant rate on every wall (not in cells). Used for the β background plasmodesmata (PD) production term in Eq. 2.

Modeled after `Degradation::OneWall` (which already existed); changed `-=` to `+=`.

---

### `Degradation::FromType`
**File:** [src/degradation.cc](src/degradation.cc), [src/degradation.h](src/degradation.h)

Decays a molecule only in cells matching a specified type. Used for extra auxin sink decay in the base ring (type = 2).

Modeled after `Degradation::One`; added a type-check condition.

---

### `PINSaturatingTransport`
**File:** [src/transport.cc](src/transport.cc), [src/transport.h](src/transport.h)

Directional, Michaelis-Menten–saturating PIN-mediated auxin transport between neighbouring cells:

```
dA_i/dt += T * Σ_j [ Pji * Aj/(1+Aj) − Pij * Ai/(1+Ai) ]
```

Also saves the signed net flux per wall into a wall variable pair (flux_ij, flux_ji) so that the WTF term in `MembraneCycling::UTGWTF` can read it.

Modeled after `DiffusionActiveTransportCell`; replaced linear transport with Michaelis-Menten saturation and added flux-save logic.

---

### `MembraneCycling::UTGWTF`
**File:** [src/membraneCycling.cc](src/membraneCycling.cc), [src/membraneCycling.h](src/membraneCycling.h)

Combined Up-The-Gradient (UTG) + With-The-Flux (WTF) PIN membrane allocation (Eq. 4):

```
dPij/dt = k_U * Pi * f(Aj) / (1 + Pi)            [UTG]
        + Pi / (1 + Pi) * (k_Wq * flux² + k_Wl * flux)  [WTF]
        − k_off * Pij                              [recycling]
```

where `f(x) = K*x / (K + x)` (saturating Hill function) and `flux` is the saved PIN-transport flux read from the wall variable written by `PINSaturatingTransport`.

Modeled after `MembraneCycling::CellUpTheGradientNonLinear`; added the WTF term (reads saved flux) and changed the recycling from a separate endo parameter to `k_off * Pij`.

---

### Factory registrations
**File:** [src/baseReaction.cc](src/baseReaction.cc)

Four `else if` blocks added to the reaction factory function to map the string names used in `.model` files to the new classes:
- `"Creation::FromType"` → `Creation::FromType`
- `"Creation::OneWall"` → `Creation::OneWall`
- `"Degradation::FromType"` → `Degradation::FromType`
- `"PINSaturatingTransport"` → `PINSaturatingTransport`
- `"MembraneCycling::UTGWTF"` → `MembraneCycling::UTGWTF`

---

## 2. Tissue Generator — `generate_tissue_files2.py`

**File:** [examples/tutorials/semicircle_leaf2/generate_tissue_files2.py](examples/tutorials/semicircle_leaf2/generate_tissue_files2.py)

This script generates all three simulation input files (`.init`, `.model`, `.solver`) each time it is run. The `.model` and `.solver` files should be edited here, not directly, because re-running the script overwrites them.

### What changed from the original:

**Tissue geometry:**
- Changed from `nRings=4, nCols=14` (56 cells) to `nRings=6, nCols=20` (120 cells)
- Added flat base-row vertices to prevent `MoveVertexRadially` from distorting the flat tissue edge
- Fixed semicircle vertex layout so inner cells form true arcs instead of squashed shapes

**Cell variables (completely renumbered):**

| Index | Old | New |
|---|---|---|
| 0 | volume tracker | volume tracker (same) |
| 1 | type | type (same) |
| 2 | (various legacy) | Aprec — auxin precursor |
| 3 | auxin A | auxin A (same) |
| 4 | PIN cytoplasmic P | PIN cytoplasmic P (same) |

**Wall variables (renumbered and expanded):**

| Index | Old | New |
|---|---|---|
| 0 | length (built-in) | length (built-in) |
| 1, 2 | Pij, Pji | Pij, Pji (same) |
| 3 | Dij | Dij (PD area) |
| 4 | (unused) | Dij mirror slot |
| 5, 6 | (absent) | flux_ij, flux_ji — saved transport flux |

**Model reactions (`.model` string):** Completely rewritten to implement the full Holloway et al. 2025 Eqs. 1–4. Old model used `DiffusionActiveTransportCell` + `AuxinModel1S`-style single-step combined rule. New model separates every biological term into its own reaction class so parameters are transparent.

**Solver:** End time reduced from `t=200` to `t=112` to stop before the second mass-division wave (which causes the RK5 adaptive solver to time out on a 160+ cell tissue).

---

## 3. New and Modified Tutorial Files

### `examples/tutorials/semicircle_leaf2/model_guide.html` — NEW
An HTML reference guide explaining the model, how to run the simulation, and what each ParaView variable represents. Includes the correct run command argument order.

### `examples/tutorials/semicircle_leaf2/semicircle_leaf2.model.annotated.txt` — NEW
A line-by-line annotated copy of the current `.model` file. Lists every reaction's source file, parameter names, parameter values, and which cell/wall variable index each number refers to. (The actual `.model` file cannot contain comments; this companion file is for reference only.)

### `examples/tutorials/semicircle_leaf2/semicircle_leaf2.model` — MODIFIED
Regenerated by `generate_tissue_files2.py`. Now contains 14 reactions + 2 compartment changes implementing the full PIN+PD model, up from a simpler 6-reaction model.

### `examples/tutorials/semicircle_leaf2/semicircle_leaf2.solver` — MODIFIED
End time reduced from 200 → 112; checkpoint interval reduced from 100 → 56 (one checkpoint per 2 time units).

---

## 4. Dij Symmetry Fix

Plasmodesmata pores are physical channels through the cell wall — the same pore seen from both sides, so Dij must equal Dji at all times. The original implementation tracked Dij in wall slot 3 but left the mirror slot (wall slot 4) static at its initial value of 0.5.

### C++ changes

**`src/transport.cc` — `DiffusionConductiveSimple::derivs`**
- Constructor relaxed: previously required exactly 1 wall index at level 1; now accepts 1 or 2
- `derivs()`: computes the conductance derivative once, writes it to the primary slot, and — when a second index is provided — writes the identical value to the mirror slot

**`src/creation.cc` — `Creation::OneWall::derivs`**
- Constructor relaxed: previously required exactly 1 wall index at level 0; now accepts 1 or 2
- `derivs()`: writes the background production rate `beta_pd` to both wall slots when a second index is present

Both changes are backwards-compatible: existing `.model` files with a single wall index continue to work unchanged.

### Model file change (`generate_tissue_files2.py`)

- `Diffusion::ConductiveSimple` index spec: `5 2 1 1` → `5 2 1 2` (second index = wall[4] = Dij_mirror)
- `Creation::OneWall` index spec: `1 1 1` → `1 1 2` (second index = wall[4] = Dij_mirror)

After this fix, wall variable 1 in ParaView shows two equal components instead of one dynamic and one frozen.

---

## 5. Saved Simulation Results

### `examples/tutorials/semicircle_leaf2/results/full_pinpd_model_t100_4x14/`
Output from an early 4×14 tissue run to t=100. Useful as a before-division baseline.

### `examples/tutorials/semicircle_leaf2/results/full_pinpd_model_t112_6x20/`
Output from a 6×20 run to t=112 before the Dij symmetry fix. Wall variable 1 components are asymmetric (slot 3 dynamic, slot 4 frozen at 0.5).

### `examples/tutorials/semicircle_leaf2/results/full_pinpd_t112_6x20_dij_symmetric/`
**Current reference run.** 6×20 tissue, t=0 to t=112, 57 timesteps, 160 cells at end. Dij and Dij_mirror are now equal throughout. Load `tissue.pvd` in ParaView (use Extract Block or load `VTK_cells*` and `VTK_walls*` as separate file series to color cells and walls independently).

---

## 6. What was NOT changed

- All original VirtualLeaf C++ files other than `transport.cc`, `creation.cc`, `degradation.cc`, `membraneCycling.cc`, and `baseReaction.cc`
- The build system (`Makefile`) — the binary was rebuilt manually. The correct WSL link command is:
  ```
  g++ -O3 -o ../bin/simulator $(ls *.o) simulator/simulator.o ply/ply_parser.o -lpthread
  ```
  `ply/ply_parser.o` must be included explicitly — it provides `ply::ply_parser::parse()` which `ply_reader.o` depends on but Boost (required to rebuild it) is not installed in WSL.
- All other tutorial examples in `examples/tutorials/`

---

## Quick reference: parameter values

| Parameter | Symbol in paper | Value used | Paper value |
|---|---|---|---|
| Radial growth rate | — | 0.006 | N/A |
| PIN transport rate | T | 3.0 | 6.0 |
| Auxin production | auxpr | 0.05 | ~0.1 |
| Auxin decay | auxdec | 0.01 | ~0.01 |
| Sink extra decay | — | 0.15 | high |
| PD diffusion | D | 0.15 | 0.8 |
| PD flux feedback | α | 0.02 | — |
| PD background prod. | β | 0.02 | — |
| PD decay ratio | γ/α | 5.0 | — |
| PIN synthesis | pinpr | 0.05 | ~0.05 |
| PIN decay | pindec | 0.02 | ~0.02 |
| UTG rate | k_U | 0.1 | 4×10⁻³ (rescaled) |
| WTF quadratic | k_Wq | 0.15 | ~0.5 (rescaled) |
| WTF linear | k_Wl | 0.2 | 3×10⁻³ (rescaled) |
| PIN recycling | k_off | 0.05 | ~0.05 |
| Division threshold | V_threshold | 6.0 | N/A |
