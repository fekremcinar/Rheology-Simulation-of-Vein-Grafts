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
| `ū` | Cross-sectional mean (average) axial velocity | m/s |
| `U_center(t)` | Centreline (peak) axial velocity at time _t_ | m/s |
| `U_mean` | Cross-sectional mean velocity (same as ū) | m/s |
| `U_max` | Peak systolic centreline velocity (~0.5 m/s) | m/s |
| `U_min` | End-diastolic centreline velocity (~0.05 m/s) | m/s |
| `v(r)` | Axial velocity as a function of radial position | m/s |
| `v_max` | Centreline (maximum) velocity in Poiseuille flow; `v_max = 2ū` | m/s |
| `p` | Kinematic pressure field (`p/ρ` in icoFoam) | m²/s² |
| `ΔP` | Pressure drop between inlet and outlet (`P_in − P_out`) | Pa |
| `Q` | Volumetric flow rate (`Q = ū · πR²`) | m³/s |
| `ρ` | Blood density (1060 kg/m³) | kg/m³ |
| `μ` | Dynamic viscosity of blood (0.0035 Pa·s at 37°C) | Pa·s |
| `ν` | Kinematic viscosity (`ν = μ/ρ ≈ 3.3×10⁻⁶`) | m²/s |
| `Re` | Reynolds number (`Re = ū·D/ν`) | — |
| `T` | Cardiac cycle period (0.857 s at 70 bpm) | s |
| `t` | Simulation time | s |
| `r` | Radial distance from vessel axis (`r = √(y²+z²)`) | m |
| `WSS` | Wall Shear Stress (`WSS = 4μū/R` for Poiseuille flow) | Pa |
| `L_entry` | Hydrodynamic entry length (`≈ 0.06·Re·D`) | m |
| `τ` | Viscous diffusion timescale (`τ = R²/ν`) | s |
| `E` | Young's modulus of vessel wall material | Pa |
| `ν_s` | Poisson's ratio of solid wall (≈ 0.45) | — |
| `t_wall` | Vessel wall thickness | m |
| `D` | Vessel diameter (`D = 2R`) | m |
| `CFL` | Courant–Friedrichs–Lewy number (`= U·deltaT/deltaX`) | — |
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

**Mean velocity (ū)**
: The cross-sectional average of the axial velocity, equal to the volumetric flow rate divided by the cross-sectional area: `ū = Q / (πR²)`. For Poiseuille flow, `ū = v_max / 2`.

**Pressure drop (ΔP)**
: The difference in static pressure between the inlet and outlet of a flow segment: `ΔP = P_in − P_out`. For fully-developed Poiseuille flow: `ΔP = 8μūL/R²`. Note that OpenFOAM's `icoFoam` stores kinematic pressure `p = P/ρ`, so the actual pressure drop in Pa is `ΔP = (p_in − p_out) × ρ`.

**Volumetric flow rate (Q)**
: Volume of fluid passing through a cross-section per unit time: `Q = ū · πR²`. For Poiseuille flow: `Q = πR⁴ΔP / (8μL)`. Units: m³/s (or mL/s in physiological contexts).

**Viscous diffusion timescale (τ)**
: The characteristic time for viscous effects to propagate radially across the vessel: `τ = R²/ν`. In a transient simulation starting from rest, the flow approaches the fully-developed Poiseuille profile after roughly `t > τ`. For Experiment 01: `τ = (0.005)²/3.3×10⁻⁶ ≈ 7.6 s`.

**Hydrodynamic entry length (L_entry)**
: The axial distance from the inlet required for the velocity profile to transition from the inlet condition (e.g., uniform plug flow) to the fully-developed Poiseuille parabola. For laminar pipe flow: `L_entry ≈ 0.06 · Re · D`. If the tube is shorter than `L_entry`, the profile is still developing at the outlet.

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

**Finite Volume Method (FVM)**
: The numerical technique used by OpenFOAM. The domain is divided into small polyhedral control volumes (cells). Conservation equations (mass, momentum) are integrated over each cell and converted to algebraic equations solved iteratively.

**Control volume (cell)**
: A small, fixed region of space in the mesh. Each cell stores field values (velocity, pressure) at its geometric centre. Fluxes are computed across cell faces.

**Face**
: A flat polygon forming the boundary between two adjacent cells (internal face), or between a cell and a boundary patch (boundary face). Used in OpenFOAM to compute fluxes between cells.

**Face centre**
: The geometric centroid of a mesh face. Boundary field values (e.g., velocity at the outlet) are associated with face centres of boundary faces, giving one sample point per boundary cell.

**Hexahedral cell (hex cell)**
: A cell with 6 quadrilateral faces and 8 vertices — the 3D equivalent of a cube. Preferred in structured meshes for accuracy and convergence. All cells in Experiments 01–07 are hexahedral.

**Mesh grading**
: Non-uniform spacing of cells within a block, controlled by an expansion ratio. Radial grading 4:1 toward the wall (used in Experiments 01–07) places smaller cells near the wall to resolve the steep velocity gradient in the boundary layer.

**Structured mesh**
: A mesh where cells are arranged in a regular, logically Cartesian grid (i, j, k indices). Offers better numerical accuracy than unstructured meshes for simple geometries. All experiments in this project use structured meshes generated by `blockMesh`.

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

**Boundary condition (BC)**
: A mathematical constraint applied at the boundary of the simulation domain (inlet, outlet, wall) that tells the solver the value or behaviour of a field at that surface. Common types used in this project: `fixedValue` (prescribes an exact value), `zeroGradient` (zero normal derivative — lets the field float freely), `noSlip` (sets velocity to zero at a wall), and `codedFixedValue` (user-defined expression compiled at runtime).

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

**PISO / SIMPLE / PIMPLE algorithms**
: Iterative pressure-velocity coupling algorithms used to solve the incompressible Navier-Stokes equations:
- **SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations) — steady-state solver.
- **PISO** (Pressure-Implicit Split-Operator) — transient solver; used inside `icoFoam`.
- **PIMPLE** — hybrid of PISO and SIMPLE; allows larger time steps with outer iterations; used inside `pimpleFoam`.

**StreamTracer**
: A ParaView filter that traces streamlines through a vector field (e.g., velocity `U`). Reveals flow patterns, recirculation zones, and secondary flows. Primary visualisation tool across all experiments.

**Glyph**
: A ParaView filter that renders arrows (or other shapes) at each point, oriented and scaled by a vector field. Used to display velocity direction and magnitude on a cross-section.

**Plot Over Line**
: A ParaView filter that samples a field along a straight line between two points and displays the result as an XY chart. Used to extract and verify the radial velocity profile (Hagen-Poiseuille parabola).

**Clip**
: A ParaView filter that cuts the geometry with a plane and displays only one half. Used to look inside the cylindrical vessel.

**Slice**
: A ParaView filter that extracts a 2D cross-section through the 3D domain along a specified plane. Used to visualise the velocity colour map on a mid-plane.

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

