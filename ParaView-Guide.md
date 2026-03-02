# Visualising Results in ParaView

This guide covers both the **one-click macro** (recommended) and the **manual step-by-step** workflow for exploring OpenFOAM simulation results in ParaView.

Each experiment ships with a dedicated Python macro under `assets/paraview/` that builds the full pipeline automatically. The manual steps below are provided for reference and for cases where you want to explore filters interactively.

---

## One-click Macro (Recommended)

Every experiment has a pre-built macro that opens the correct case, applies all filters, and arranges the layout automatically.

**How to register a macro (one-time setup per experiment):**

1. Open **ParaView** on macOS.
2. **Tools → Macros → Add new macro** → select the script for your experiment → **OK**.
3. The macro then appears in the **Macros** menu — click it whenever you want to reload the visualisation.

| Experiment | Macro script |
|---|---|
| 01 — Simple Laminar Flow | `assets/paraview/01_simple_laminar.py` |
| 02 — Heartbeat Laminar Flow | `assets/paraview/02_heartbeat_laminar.py` |

> **Note:** Each time the macro runs it first closes **all currently open views and pipeline objects without saving**, then rebuilds everything from scratch. This ensures a clean session on every run.

---

## Manual Step-by-Step (Experiment 01 — Simple Laminar Flow)

The steps below are specific to `01_simple_laminar`. The same filters and logic apply to every experiment; adjust the field names and coordinates accordingly.

### Step 1 — Open the case

1. On macOS (outside Docker), open ParaView.
2. **File → Open** → navigate to `~/Rheology-Simulation-of-Vein-Grafts/run/01_simple_laminar/`.
3. Select `01_simple_laminar.foam` → **OK** → click **Apply**.
   You will see a full 3D cylinder in the viewport.

### Step 2 — Go to the last time step

1. Click the **⏭ Last Frame** button in the toolbar to jump to t = 10 s (fully developed flow).

### Step 3 — Colour by axial velocity

1. Change the field dropdown from `p` to **U**, component **X**.

You will see **blue at the wall** (U = 0, no-slip) grading to **red at the centre** (peak ≈ 0.2 m/s).

### Step 4 — Clip the cylinder to see inside

1. Select `01_simple_laminar.foam` in the Pipeline Browser.
2. **Filters → Search** (or press **Space**) → type `Clip` → select **Clip** 
3. Set **Normal = (0, 0, 1)**, **Origin = (0, 0, 0)** to remove the top half.
4. Click **Apply**, you will see the inside of the cylinder coloured by velocity.

### Step 5 — Slice down the middle to see the full profile

1. Select `01_simple_laminar.foam` in the Pipeline Browser.
2. **Filters → Search** → type `Slice` → select **Slice**
3. Set **Normal = (0, 0, 1)**, **Origin = (0, 0, 0)**.
   Colour by **U → X**.
4. Click **Apply**. You will see the parabolic colour gradient: blue at both walls (U = 0), red at the centreline (U ≈ 0.2 m/s).

> **Note:** The parabola is shown as a colour gradient on the slice surface, not as a drawn curve.

### Step 6 — Verify the parabolic profile numerically (Plot Over Line)

1. Select `01_simple_laminar.foam` in the Pipeline Browser (important: select the source, not Clip or Slice).
2. **Filters → Search** → type `Plot Over Line` → select it.
3. Set:
   - **Point 1**: `(0.09, -0.005, 0)` — near outlet, bottom wall
   - **Point 2**: `(0.09,  0.005, 0)` — near outlet, top wall
4. Click **Apply**. A **LineChartView** opens on the right. In the chart series list, **uncheck everything except `U_X`**.
   You should see a smooth parabola: zero at both walls (y = ±0.005), peak ≈ 0.2 m/s at the centre (y = 0).

### Step 7 — Streamlines (standard method for all experiments)

Streamlines are the primary visualization tool across all experiments. They naturally reveal recirculation zones, flow separation, and secondary flows in later experiments (junctions, grafts).

1. Select `01_simple_laminar.foam` in the Pipeline Browser.
2. **Filters → Search** (or press **Space**) → type `Stream Tracer` → select it.
3. Set in Properties:
   - **Vectors**: `U`
   - **Seed Type**: `Line Source`
   - **Point 1**: `(0.001, -0.004, 0)` — near inlet, close to wall
   - **Point 2**: `(0.001,  0.004, 0)` — near inlet, opposite side
   - **Resolution**: `20`
   - **Max Streamline Length**: `0.2`
4. colour by **U → Magnitude**. Click **Apply**

In a straight tube the streamlines will be straight and parallel. In future junction and graft experiments the same setup will automatically reveal recirculation zones and spiral flows.

### Step 8 — Velocity arrows on a cross-section (Glyph)

1. Select `Slice1` in the Pipeline Browser.
2. **Filters → Search** → type `Glyph` → select it.
3. Set:
   - **Glyph Type**: `Arrow`
   - **Orientation Array**: `U`
   - **Scale Array**: `U`
   - **Scale Factor**: `0.02`
   - **Glyph Mode**: `Every Nth Point`, N = `2`
4. Click **Apply**.
   Arrows show longer in the centre (fast) and shorter/zero at the walls.

---

## General Tips

| Task | Filter / Tool |
|---|---|
| Expose vessel interior | **Clip** — Normal (0,0,1), Origin (0,0,0) |
| Longitudinal colour profile | **Slice** — Normal (0,0,1), colour by U_X |
| Flow paths / recirculation | **StreamTracer** — seed line near inlet |
| Velocity magnitude arrows | **Glyph** — Arrow type, scale by U |
| Radial velocity profile chart | **Plot Over Line** — vertical line at chosen x |
| Derived quantities (WSS, etc.) | **Calculator** filter |

ParaView User's Guide: <https://docs.paraview.org/en/latest/UsersGuide/index.html>
