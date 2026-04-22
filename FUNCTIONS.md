# Tissue Simulator — Reaction Rules Reference

All reaction rules used in this project, what they do, and their parameters.

---

## AuxinModel1S

The main rule. Handles auxin transport and PIN protein cycling together in one combined rule. This is the same rule used in the Caggiano & Jonsson (2017) KAN/REV paper.

**What it does biologically:**
- Auxin is produced in all cells at a low constant rate and decays over time
- Auxin moves between neighbouring cells two ways: slowly by passive diffusion, and actively by PIN protein pumps
- PIN proteins are produced in proportion to local auxin levels
- PIN preferentially accumulates on the wall face pointing toward whichever neighbour has more auxin (up-the-gradient polarization)
- This positive feedback between auxin flow and PIN positioning causes auxin to self-organize into a narrow stripe (the future vein)

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| c_A | 0.005 | How fast auxin is produced in all cells |
| d_A | 0.005 | How fast auxin decays |
| k_p | 0.9 | Fraction of PIN that is polarized — 0.9 means 90% directional, 10% symmetric |
| f_p | 0.3 | Strength of PIN polarization feedback |
| T | 1.3 | Active transport rate through PIN — higher = stronger canalization |
| D | 0.002 | Passive diffusion rate — kept low so PIN dominates transport |
| c_P | 0.005 | How fast PIN is produced |
| d_P | 0.005 | How fast PIN decays |

**Cell variables used:**

| Variable | What it stores |
|---|---|
| cell variable 4 | Auxin concentration |
| cell variable 5 | PIN cytoplasm pool |
| cell variable 6 | X boundary marker (not used in this model) |
| cell variable 7 | Membrane marker |

**Wall variable used:**

| Variable | What it stores |
|---|---|
| wall variable 1 | PIN on membrane (one value per wall side) |

---

## WallMechanics::Spring

Gives the tissue its mechanical structure by making cell walls behave like springs.

**What it does biologically:**
Each cell wall has a natural resting length. If a wall gets stretched or compressed, it pushes back to restore that length. This keeps cells from collapsing or exploding and maintains the overall tissue shape.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| k_force | 0.01 | Spring stiffness — how strongly walls resist deformation |
| frac_adhesion | 1.0 | Extra resistance when walls are compressed (cell adhesion) |

---

## CenterCOM

Keeps the tissue centred on screen.

**What it does:**
At every time step, calculates the centre of mass of all cells and moves everything so the centre sits at coordinates (0, 0). Without this the tissue would slowly drift off screen due to uneven forces.

**Parameters:** None

---

## MoveVertexRadially

Makes the tissue grow by pushing all cell vertices outward from the centre.

**What it does biologically:**
Simulates turgor-driven radial growth. Each vertex is pushed outward at a rate proportional to its distance from the centre, producing exponential growth where outer cells grow faster than inner cells.

**Note:** This rule pushes ALL vertices including the flat base edge, which can cause the base corners to distort over long simulation times. This is a known cosmetic limitation.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| k_growth | 0.001 | Growth rate — how fast the tissue expands |
| r_pow | 1 | 1 = exponential growth, 0 = linear growth |

---

## WallGrowth::Stress

Makes walls grow permanently when they are stretched, simulating plant cell wall yielding.

**What it does biologically:**
When a wall is stretched beyond its resting length by mechanical stress, the resting length is permanently increased — the wall has grown. This is based on the Lockhart model of plant growth and is more biologically realistic than simply pushing vertices outward.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| k_growth | 0.005 | How fast walls grow under stress |
| stress_threshold | 0.0 | Minimum stress needed before growth starts |
| stretch_flag | 1 | 1 = respond to strain (relative stretch), 0 = respond to stress variable |
| linear_flag | 1 | 1 = growth proportional to wall length |

---

## Creation::Zero

Produces a molecule at a constant rate in every cell.

**What it does:**
Adds a fixed amount of a molecule to every cell at every time step, regardless of what else is happening. Used for background auxin production.

**Parameters used:**

| Parameter | What it controls |
|---|---|
| k_production | How much is produced per time step |

---

## Degradation::One

Breaks down a molecule at a rate proportional to how much is present (first-order decay).

**What it does:**
Removes a fraction of the molecule every time step. Molecules at high concentration decay faster than at low concentration, eventually reaching a steady state when production equals decay.

**Parameters used:**

| Parameter | What it controls |
|---|---|
| k_decay | Fraction removed per time step |

---

## DiffusionSimple

Moves a molecule between neighbouring cells by passive diffusion.

**What it does biologically:**
Molecules flow from high concentration cells to low concentration neighbours, driven purely by concentration difference. Does not cross the boundary edge of the tissue. Used for testing before adding PIN-mediated transport.

**Parameters used:**

| Parameter | What it controls |
|---|---|
| diffusion_rate | How fast molecules move between cells |

---

## DiffusionActiveTransportCell

Combines passive diffusion with PIN-mediated active transport between cells.

**What it does:**
Two transport mechanisms work together:
- Passive diffusion (rate D) — moves auxin from high to low concentration
- Active transport (rate T) — PIN pumps on wall faces actively move auxin in a specific direction, which can work against the concentration gradient

Setting T=0 gives pure diffusion. Setting T high gives strong PIN-mediated canalization.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| D | 0.1 | Passive diffusion rate |
| T | 2.0 | Active PIN transport rate |

---

## MembraneCycling::CellUpTheGradientNonLinear

Controls how PIN proteins move between the cell interior and the cell membrane (up-the-gradient version).

**What it does biologically:**
PIN is constantly cycling between the cytoplasm (inside the cell) and the membrane (on the cell wall). This rule makes PIN preferentially insert into the membrane face pointing toward whichever neighbour has more auxin. The response follows a Hill function so even small differences in neighbour auxin create strong PIN polarization.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| k_exo | 0.1 | Rate of PIN insertion into membrane |
| k_endo | 0.05 | Rate of PIN removal from membrane |
| K | 1.0 | Auxin concentration needed for half-maximal PIN insertion |
| n | 2.0 | Steepness of response — higher = more switch-like |

---

## Division::ShortestPath2D

Divides a cell when it reaches a size threshold, cutting along the shortest possible path.

**What it does biologically:**
When a cell grows to twice its original size, it divides. The division plane is chosen to minimize the length of the new wall created, following the principle of least mechanical work. Only works in 2D tissue geometries.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| V_threshold | 500.0 | Cell size at which division is triggered |
| LWall_frac | 1.0 | Length of new wall relative to actual cut length |
| Lwall_threshold | 0.3 | Minimum distance from existing vertices (avoids degenerate geometry) |
| COM_flag | 1 | 1 = divide through centre of mass |

---

## RemovalOutsideRadius

Removes cells that grow beyond a specified radius from the origin.

**What it does:**
Acts as a boundary condition — any cell whose centre moves beyond the radius threshold is deleted from the simulation. Used to keep the tissue within a defined region during growth.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| R_threshold | 55.0 | Radius beyond which cells are removed |

---

## References

- Jonsson H, Heisler M, Shapiro B, Meyerowitz E, Mjolsness E (2006). An auxin-driven polarized transport model for phyllotaxis. *PNAS* 103:1633–1638.
- Holloway DM, Eiriksson TK, Wenzel CL (2025). The role of auxin transport through plasmodesmata in leaf vein canalization and patterning. *Frontiers in Plant Science* 16:1621815.
- Caggiano MP et al. (2017). Cell type boundaries organize plant development. *eLife* 6:e27421.
- Merks RMH et al. (2011). VirtualLeaf. *Plant Physiology* 155:656–666.