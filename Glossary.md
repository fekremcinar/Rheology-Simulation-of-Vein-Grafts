# Glossary — Rheology Simulation of Vein Grafts

A reference of technical and scientific terms, symbols, and annotations used throughout this project.

---

## Symbols

| Symbol | Description | Unit |
|--------|-------------|------|
| `r1` | Donor (upstream) vessel radius | m |
| `r2` | Recipient (downstream) vessel radius | m |
| `r3` | Venous graft radius (intermediate bridge) | m |
| `R` | Generic vessel radius in formulae | m |
| `L` | Venous graft length | m |
| `U` | Velocity vector field | m/s |
| `U_center(t)` | Centreline (peak) axial velocity at time _t_ | m/s |
| `U_mean` | Cross-sectional mean velocity | m/s |
| `U_max` | Peak systolic centreline velocity (~0.5 m/s) | m/s |
| `U_min` | End-diastolic centreline velocity (~0.05 m/s) | m/s |
| `p` | Kinematic pressure field (`p/ρ` in icoFoam) | m²/s² |
| `ρ` | Blood density (1060 kg/m³) | kg/m³ |
| `μ` | Dynamic viscosity of blood (0.0035 Pa·s at 37°C) | Pa·s |
| `ν` | Kinematic viscosity (`ν = μ/ρ ≈ 3.3×10⁻⁶`) | m²/s |
| `Re` | Reynolds number (`Re = U·D/ν`) | — |
| `T` | Cardiac cycle period (0.857 s at 70 bpm) | s |
| `t` | Simulation time | s |
| `r` | Radial distance from vessel axis (`r = √(y²+z²)`) | m |
| `WSS` | Wall Shear Stress | Pa |
| `E` | Young's modulus of vessel wall material | Pa |
| `ν_s` | Poisson's ratio of solid wall (≈ 0.45) | — |
| `t_wall` | Vessel wall thickness | m |
| `D` | Vessel diameter (`D = 2R`) | m |
| `CFL` | Courant–Friedrichs–Lewy number | — |
| `deltaT` | Simulation time step | s |

---

## Medical & Physiological Terms

**Anastomosis**
: Surgical connection between two vessels, ducts, or other tubular structures. In this project it refers to the end-to-end suture of the donor and recipient vessel ends.

**Venous graft**
: A segment of vein (or synthetic equivalent) used as an intermediate bridging conduit between two vessels of significantly different radii. Required when `r1/r2 > 3/2` or `r1/r2 < 2/3` to prevent severe flow disturbance.

**Donor vessel (artery)**
: The upstream vessel supplying blood to the junction (radius `r1`).

**Recipient vessel (artery)**
: The downstream vessel receiving blood from the junction (radius `r2`).

**Systole**
: The phase of the cardiac cycle in which the heart contracts and ejects blood. Corresponds to peak velocity (`U_max`).

**Diastole**
: The phase of the cardiac cycle in which the heart relaxes and refills. Corresponds to minimum velocity (`U_min`).

**Heart rate**
: Number of cardiac cycles per minute. Baseline value used in this project: 70 bpm (period T = 0.857 s).

**Pulsatile flow**
: Flow driven by a periodic pressure gradient mimicking the heartbeat, as opposed to steady (constant) flow. Implemented in Experiments 02–07 via a sinusoidal centreline velocity waveform.

**Wall Shear Stress (WSS)**
: Tangential stress exerted by the flowing blood on the vessel wall. Physiological range for healthy endothelium: **0.5–4.0 Pa**. Both too-low WSS (favours atherosclerosis) and too-high WSS (endothelial damage) are pathological.

**Endothelium**
: The thin layer of cells lining the inner surface of blood vessels. Highly sensitive to WSS; abnormal WSS promotes vascular disease.

**Atherosclerosis**
: Accumulation of plaques in artery walls, promoted by chronically low or oscillatory WSS (< 0.5 Pa).

**Compliance (vascular)**
: The ability of a vessel wall to expand and recoil in response to pulsatile pressure changes. Venous walls are more compliant (lower _E_) than arterial walls.

---

## Fluid Mechanics Terms

**Laminar flow**
: Smooth, orderly flow in parallel layers with no lateral mixing. Characterised by Reynolds number `Re < 2300` in pipe flow.

**Turbulent flow**
: Chaotic, irregular flow with lateral mixing. Begins at `Re > ~4000` in pipe flow. Pathological in vascular grafts.

**Reynolds number (Re)**
: Dimensionless ratio of inertial to viscous forces: `Re = U·D/ν`. Determines whether flow is laminar or turbulent.

**Poiseuille flow (Hagen-Poiseuille)**
: Exact analytical solution for fully-developed, steady, laminar flow in a straight circular tube. Velocity profile is a parabola: `U(r) = U_center · (1 − r²/R²)`. Peak centreline velocity is twice the mean velocity.

**Parabolic velocity profile**
: The cross-sectional velocity distribution in Poiseuille flow — zero at the wall (no-slip) and maximum at the centreline.

**No-slip condition**
: Boundary condition requiring fluid velocity to match the wall velocity. For a stationary wall: `U_wall = 0`.

**Recirculation zone**
: A region of reversed or separated flow, typically occurring downstream of an abrupt change in geometry (e.g., a step junction). Associated with low or negative WSS. The primary hazard in mismatched vessel junctions.

**Reattachment length**
: The axial distance downstream of a step at which the separated flow reattaches to the wall. Used in Experiment 06 to evaluate graft length optimality.

**Entry length**
: The axial distance required for a developing velocity profile to reach the fully developed (parabolic) state. For laminar pipe flow: `L_entry ≈ 0.06 · Re · D`.

**Courant–Friedrichs–Lewy (CFL) number**
: Numerical stability criterion for explicit time-stepping: `CFL = U·deltaT/deltaX`. Must be `< 1` for stability in `icoFoam`.

**Fluid-Structure Interaction (FSI)**
: Coupled simulation of fluid flow and structural deformation of the vessel wall. Used in Experiment 03.

**Incompressible flow**
: Flow in which density is assumed constant. Blood is treated as incompressible in all experiments (`icoFoam` / `pimpleFoam`).

**Kinematic pressure**
: Pressure divided by density (`p/ρ`), as used internally by OpenFOAM's `icoFoam` solver. Units: m²/s².

**Viscosity**
: Resistance of a fluid to deformation/shear.
- *Dynamic viscosity* (`μ`): force per unit area per unit velocity gradient; units Pa·s.
- *Kinematic viscosity* (`ν = μ/ρ`): dynamic viscosity normalised by density; units m²/s.

**Newtonian fluid**
: A fluid in which viscosity is independent of shear rate. Blood is approximated as Newtonian in this project (valid for large vessels, `Re ≫ 1`).

**Rheology**
: The study of the flow and deformation of matter, particularly non-Newtonian fluids. In the biomedical context, it encompasses blood flow behaviour, including viscosity, elasticity, and the effects of flow geometry on vessel wall interactions.

---

## OpenFOAM Terms

**blockMesh**
: OpenFOAM utility that generates a structured hexahedral mesh from a `blockMeshDict` file.

**checkMesh**
: OpenFOAM utility that reports mesh quality metrics (orthogonality, aspect ratio, skewness).

**icoFoam**
: OpenFOAM solver for incompressible, laminar, transient (unsteady) Newtonian flow. Used in Experiments 01–02, 04–07.

**pimpleFoam**
: OpenFOAM solver for incompressible, laminar or turbulent, transient flow using the PIMPLE (merged PISO–SIMPLE) algorithm. Used in Experiment 03 (FSI).

**blockMeshDict**
: OpenFOAM dictionary file specifying mesh vertices, blocks, edges, boundary patches, and grading.

**controlDict**
: OpenFOAM dictionary that controls simulation time (`startTime`, `endTime`, `deltaT`), output frequency (`writeInterval`), and solver selection (`application`).

**fvSchemes**
: OpenFOAM dictionary specifying the discretisation schemes for gradient, divergence, and Laplacian operators.

**fvSolution**
: OpenFOAM dictionary specifying linear solver settings (solver type, tolerances, relaxation factors) for each field.

**transportProperties**
: OpenFOAM dictionary in `constant/` that defines fluid properties. For a Newtonian fluid: kinematic viscosity `nu`.

**codedFixedValue**
: OpenFOAM boundary condition that allows arbitrary C++ expressions to be compiled at runtime and applied as inlet/outlet values. Used for the pulsatile parabolic inlet profile in Experiments 02–07.

**O-grid (butterfly mesh)**
: A structured mesh topology for circular cross-sections: a central square (or rectangular) block surrounded by four outer curved blocks. Avoids degenerate cells at the axis that appear in purely polar meshes.

**Wedge mesh**
: A 2D-axisymmetric mesh representation in OpenFOAM using a single-cell-thick wedge (typically 5° sector) with `wedge` boundary conditions on the front and back faces.

**polyMesh**
: The directory inside `constant/` where OpenFOAM stores the generated mesh files (`points`, `faces`, `cells`, `boundary`).

**ParaView**
: Open-source, multi-platform data analysis and visualisation application. Used to post-process and visualise OpenFOAM simulation results.

---

## Clinical Context

**Tissue transplantation**
: Surgical transfer of tissue (skin, bone, organ) from one site to another, often requiring re-anastomosis of blood vessels of different calibres.

**Vessel repair surgery**
: Surgical reconstruction of damaged or diseased blood vessels, which may require interposing a graft when the size mismatch is too large for direct suture.

**Graft patency**
: The long-term openness and functionality of a vascular graft. Disturbed flow (recirculation, low WSS) is a primary cause of graft failure through intimal hyperplasia.

**Intimal hyperplasia**
: Thickening of the inner vessel wall layer (intima) in response to flow disturbance or surgical injury. The leading cause of late graft failure.

