# Tissue Simulator — Reaction Rules Reference

All reaction rules used in this project, what they do, and their parameters.

---

## AuxinModel1S

The main rule. Handles auxin transport and PIN protein cycling together in one combined rule. This is the same rule used in the Caggiano & Jonsson (2017) KAN/REV paper.

**Source:** [src/network.h:288](src/network.h#L288) · [src/network.cc](src/network.cc) · factory: [src/baseReaction.cc:448](src/baseReaction.cc#L448)

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

**Source:** [src/mechanicalSpring.h:70](src/mechanicalSpring.h#L70) · [src/mechanicalSpring.cc](src/mechanicalSpring.cc) · factory: [src/baseReaction.cc:99](src/baseReaction.cc#L99)

**Note:** Also accepted under the older name `VertexFromWallSpring`.

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

**Source:** [src/adhocReaction.h:562](src/adhocReaction.h#L562) · [src/adhocReaction.cc](src/adhocReaction.cc) · factory: [src/baseReaction.cc:606](src/baseReaction.cc#L606)

**What it does:**
At every time step, calculates the centre of mass of all cells and moves everything so the centre sits at coordinates (0, 0). Without this the tissue would slowly drift off screen due to uneven forces.

**Parameters:** None

---

## MoveVertexRadially

Makes the tissue grow by pushing all cell vertices outward from the centre.

**Source:** [src/growthForce.h:48](src/growthForce.h#L48) (class `GrowthForce::Radial`) · factory: [src/baseReaction.cc:258](src/baseReaction.cc#L258)

**Note:** `MoveVertexRadially` is an alias for `GrowthForce::Radial`; both names are accepted.

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

**Source:** [src/growth.h:115](src/growth.h#L115) · [src/growth.cc](src/growth.cc) · factory: [src/baseReaction.cc:51](src/baseReaction.cc#L51)

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

**Source:** [src/creation.h:31](src/creation.h#L31) · [src/creation.cc](src/creation.cc) · factory: [src/baseReaction.cc:328](src/baseReaction.cc#L328)

**What it does:**
Adds a fixed amount of a molecule to every cell at every time step, regardless of what else is happening. Used for background auxin production.

**Parameters used:**

| Parameter | What it controls |
|---|---|
| k_production | How much is produced per time step |

---

## Degradation::One

Breaks down a molecule at a rate proportional to how much is present (first-order decay).

**Source:** [src/degradation.h:35](src/degradation.h#L35) · [src/degradation.cc](src/degradation.cc) · factory: [src/baseReaction.cc:367](src/baseReaction.cc#L367)

**What it does:**
Removes a fraction of the molecule every time step. Molecules at high concentration decay faster than at low concentration, eventually reaching a steady state when production equals decay.

**Parameters used:**

| Parameter | What it controls |
|---|---|
| k_decay | Fraction removed per time step |

---

## DiffusionSimple

Moves a molecule between neighbouring cells by passive diffusion.

**Source:** [src/transport.h:81](src/transport.h#L81) · [src/transport.cc](src/transport.cc) · factory: [src/baseReaction.cc:416](src/baseReaction.cc#L416)

**What it does biologically:**
Molecules flow from high concentration cells to low concentration neighbours, driven purely by concentration difference. Does not cross the boundary edge of the tissue. Used for testing before adding PIN-mediated transport.

**Parameters used:**

| Parameter | What it controls |
|---|---|
| diffusion_rate | How fast molecules move between cells |

---

## DiffusionActiveTransportCell

Combines passive diffusion with PIN-mediated active transport between cells.

**Source:** [src/transport.h:337](src/transport.h#L337) · [src/transport.cc](src/transport.cc) · factory: [src/baseReaction.cc:428](src/baseReaction.cc#L428)

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

## DiffusionConductiveSimple

Passive diffusion between cells through a conductance variable that grows with flux (plasmodesmata dynamics).

**Source:** [src/transport.h](src/transport.h) · [src/transport.cc](src/transport.cc) · factory: [src/baseReaction.cc](src/baseReaction.cc)

**Note:** This is an existing VirtualLeaf class that was **modified** for this project to accept an optional second wall index (mirror slot).

**What it does:**
Two things at once:

1. **Auxin diffusion** (Eq. 1 PD term): `dA_i/dt += D * Dij * (Aj − Ai)` — moves auxin between neighbours proportionally to the plasmodesmata conductance Dij.
2. **PD area dynamics** (Eq. 2): `dDij/dt = α * |flux|^p2 / Dij^(p3+1) * Dij − γ/α * Dij * α` — Dij grows on walls that carry high flux (positive feedback) and decays otherwise, implementing canalization of PD channels.

Each wall is processed once (`if i < neighbour`), so Dij is naturally symmetric.

**Modified:** Now accepts 1 or 2 wall variable indices at level 1. When two indices are provided, the same conductance derivative is written to both the primary slot (Dij) and the mirror slot (Dij_mirror), keeping them equal. Backwards-compatible with single-index usage.

**Parameters:**

| Parameter | Value | What it controls |
|---|---|---|
| D | 0.15 | Baseline PD diffusion rate (multiplied by Dij) |
| α (alpha_pd) | 0.02 | Rate at which flux widens PD pores |
| p2 | 2 | Power of flux in feedback term (2 = flux-squared) |
| p3 | 0 | Power in denominator (0 = disabled) |
| γ/α (gamma_over_alpha) | 5.0 | Sets the steady-state Dij ceiling: Dij* ≈ flux^(p2/2) / sqrt(γ/α) |

---

## MembraneCycling::CellUpTheGradientNonLinear

Controls how PIN proteins move between the cell interior and the cell membrane (up-the-gradient version).

**Source:** [src/membraneCycling.h:220](src/membraneCycling.h#L220) · [src/membraneCycling.cc:429](src/membraneCycling.cc#L429) · factory: [src/baseReaction.cc:707](src/baseReaction.cc#L707)

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

**Source:** [src/compartmentDivision.h:538](src/compartmentDivision.h#L538) · [src/compartmentDivision.cc](src/compartmentDivision.cc) · factory: [src/baseCompartmentChange.cc](src/baseCompartmentChange.cc)

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

**Source:** [src/compartmentRemoval.h:44](src/compartmentRemoval.h#L44) · [src/compartmentRemoval.cc](src/compartmentRemoval.cc) · factory: [src/baseCompartmentChange.cc](src/baseCompartmentChange.cc)

**What it does:**
Acts as a boundary condition — any cell whose centre moves beyond the radius threshold is deleted from the simulation. Used to keep the tissue within a defined region during growth.

**Parameters used:**

| Parameter | Value | What it controls |
|---|---|---|
| R_threshold | 55.0 | Radius beyond which cells are removed |

---

---

## NEW CLASSES — added for Holloway et al. 2025 implementation

---

## Creation::FromType

Produces a molecule only in cells whose type variable matches a specified value.

**Source:** [src/creation.h](src/creation.h) · [src/creation.cc](src/creation.cc) · factory: [src/baseReaction.cc](src/baseReaction.cc)

**What it does:** Used to ramp up auxin precursor (Aprec) exclusively in the outer source ring. Cells with type=1 receive constant production; all other cells are unaffected.

**Parameters:**

| Parameter | Value | What it controls |
|---|---|---|
| k_c | 0.02 | Production rate in matching cells |
| type_value | 1 | Which cell type receives production |

---

## Creation::OneWall

Produces a molecule at a constant rate on every wall (not in cells).

**Source:** [src/creation.h](src/creation.h) · [src/creation.cc](src/creation.cc) · factory: [src/baseReaction.cc](src/baseReaction.cc)

**What it does:** Implements the background plasmodesmata (PD) production term β in Eq. 2. Every wall gains a small constant increment of Dij each step, preventing PD pores from closing completely on walls with no auxin flux.

**Modified:** Now accepts an optional second wall variable index at level 0. When two indices are given, the same production rate is written to both slots — used to keep the Dij mirror slot (wall[4]) equal to the primary Dij slot (wall[3]).

**Parameters:**

| Parameter | Value | What it controls |
|---|---|---|
| k_cw | 0.02 | Background production rate on all walls |

---

## Degradation::FromType

Decays a molecule at a higher rate in cells matching a specified type.

**Source:** [src/degradation.h](src/degradation.h) · [src/degradation.cc](src/degradation.cc) · factory: [src/baseReaction.cc](src/baseReaction.cc)

**What it does:** Implements the extra auxin sink at the base ring. Cells with type=2 (sink) get an additional fast decay on top of the background decay from `Degradation::One`, driving auxin to flow toward the base.

**Parameters:**

| Parameter | Value | What it controls |
|---|---|---|
| k_d | 0.15 | Extra decay rate in matching cells |
| type_value | 2 | Which cell type gets extra decay (2 = sink/base ring) |

---

## PINSaturatingTransport

Directional, Michaelis-Menten–saturating PIN-mediated auxin transport between cells.

**Source:** [src/transport.h](src/transport.h) · [src/transport.cc](src/transport.cc) · factory: [src/baseReaction.cc](src/baseReaction.cc)

**What it does biologically:**
PIN proteins on cell walls pump auxin from one cell to its neighbour. The pumping rate saturates at high auxin concentrations (Michaelis-Menten kinetics), preventing runaway accumulation. For each shared wall:

```
dA_i/dt += T * [ Pji * Aj/(1+Aj) − Pij * Ai/(1+Ai) ]
```

Also saves the directional net flux magnitude into a wall variable pair so `MembraneCycling::UTGWTF` can use it for the WTF term.

**Parameters:**

| Parameter | Value | What it controls |
|---|---|---|
| T | 3.0 | PIN transport permeability (paper value: 6.0) |

---

## MembraneCycling::UTGWTF

Combined Up-The-Gradient (UTG) and With-The-Flux (WTF) PIN membrane allocation rule (Eq. 4 of Holloway et al. 2025).

**Source:** [src/membraneCycling.h](src/membraneCycling.h) · [src/membraneCycling.cc](src/membraneCycling.cc) · factory: [src/baseReaction.cc](src/baseReaction.cc)

**What it does biologically:**
Controls how PIN cycles between the cytoplasm and each membrane face. Two complementary mechanisms:

- **UTG (Up The Gradient):** biases PIN insertion toward the neighbour with more auxin. Detects the local gradient — initiates canalization from a shallow auxin slope.
- **WTF (With The Flux):** reinforces PIN on walls that already carry net flux. Detects flow — narrows and consolidates channels once started.
- **Recycling (k_off):** removes PIN from membrane back to cytoplasm, allowing redistribution.

```
dPij/dt = k_U * Pi * f(Aj) / (1 + Pi)            [UTG]
        + Pi / (1 + Pi) * (k_Wq * flux² + k_Wl * flux)  [WTF]
        − k_off * Pij                              [recycling]
```

where `f(x) = K*x / (K + x)` and flux is the signed flux saved by `PINSaturatingTransport`.

**Parameters:**

| Parameter | Value | What it controls |
|---|---|---|
| k_U | 0.1 | UTG allocation rate — how strongly PIN follows the auxin gradient |
| k_off | 0.05 | PIN recycling rate — how fast PIN returns from membrane to cytoplasm |
| K | 2.0 | Half-saturation constant for the UTG auxin-sensing function |
| k_Wq | 0.15 | WTF quadratic flux coefficient (flux² term — sharpens channel boundaries) |
| k_Wl | 0.2 | WTF linear flux coefficient (flux¹ term — broader initial response) |

---

## References

- Jonsson H, Heisler M, Shapiro B, Meyerowitz E, Mjolsness E (2006). An auxin-driven polarized transport model for phyllotaxis. *PNAS* 103:1633–1638.
- Holloway DM, Eiriksson TK, Wenzel CL (2025). The role of auxin transport through plasmodesmata in leaf vein canalization and patterning. *Frontiers in Plant Science* 16:1621815.
- Caggiano MP et al. (2017). Cell type boundaries organize plant development. *eLife* 6:e27421.
- Merks RMH et al. (2011). VirtualLeaf. *Plant Physiology* 155:656–666.
