#!/usr/bin/env pvpython
"""
Standard ParaView visualization — Rheology Simulation of Vein Grafts
======================================================================
Applies the standard pipeline to any OpenFOAM icoFoam case in this project:
  1. Clip        — remove top half to expose the vessel interior
  2. Slice       — longitudinal midplane coloured by axial velocity (U_X)
  3. StreamTracer — streamlines (reveals recirculation in junction/graft cases)
  4. Glyph       — velocity arrows on the slice
  5. PlotOverLine — U_X profile chart near the outlet (parabola validation)

How to run
----------
Option A — pvpython (command line, outside ParaView):
    pvpython ~/Rheology-Simulation-of-Vein-Grafts/assets/paraview/standard_viz.py \\
        --case ~/Rheology-Simulation-of-Vein-Grafts/run/01_simple_laminar \\
        --radius 0.005

Option B — ParaView Python Shell (Tools → Python Shell):
    1. Edit CASE_DIR and RADIUS in the "── Macro defaults ──" section below.
    2. Open the script in the Python Shell and click Run Script.
    The pipeline panels will appear directly in the ParaView layout.

Option C — ParaView Macro (reuse with one click):
    1. Edit CASE_DIR and RADIUS in the "── Macro defaults ──" section below.
    2. Tools → Macros → Add new macro → select this file.
    3. Click the macro from the Macros menu whenever you open a new case.

Parameters
----------
--case    Path to the run case folder. Must contain <case_name>.foam.
--radius  Inner vessel radius in metres (default: 0.005).
"""

import sys
import os

# ── Detect execution context ──────────────────────────────────────────────────
# True when launched via pvpython; False when run inside the ParaView GUI
# (Tools → Python Shell or Tools → Macros).
_IS_PVPYTHON = 'pvpython' in sys.executable

# ── Macro defaults (edit these when using inside the ParaView GUI) ────────────
CASE_DIR = os.path.expanduser(
    "~/Rheology-Simulation-of-Vein-Grafts/run/01_simple_laminar"
)
RADIUS = 0.005  # vessel inner radius [m]

# ── Parse command-line arguments (overrides macro defaults; pvpython only) ────
if _IS_PVPYTHON:
    try:
        import argparse
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--case',   default=CASE_DIR,
                            help='Path to the run case folder')
        parser.add_argument('--radius', type=float, default=RADIUS,
                            help='Vessel inner radius [m]')
        args, _ = parser.parse_known_args()
        CASE_DIR = os.path.expanduser(args.case)
        RADIUS   = args.radius
    except Exception:
        pass

CASE_NAME = os.path.basename(CASE_DIR.rstrip('/'))
FOAM_FILE = os.path.join(CASE_DIR, CASE_NAME + '.foam')

# ── ParaView API ──────────────────────────────────────────────────────────────
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# ── Load the OpenFOAM case ────────────────────────────────────────────────────
print(f"\n[standard_viz] Loading: {FOAM_FILE}")
reader = OpenFOAMReader(FileName=FOAM_FILE)
reader.MeshRegions = ['internalMesh']
reader.CellArrays  = ['U', 'p']
UpdatePipeline()

# Detect mesh bounds for auto-scaling
bounds   = reader.GetDataInformation().GetBounds()
# bounds = (x_min, x_max, y_min, y_max, z_min, z_max)
x_outlet = bounds[0] + 0.9 * (bounds[1] - bounds[0])   # 90% toward outlet

print(f"[standard_viz]   Vessel radius : {RADIUS} m")
print(f"[standard_viz]   Profile line  : x = {x_outlet:.4f} m")

# Jump to last time step (fully developed flow)
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.GoToLast()
UpdatePipeline()

# ── Render view ───────────────────────────────────────────────────────────────
renderView = GetActiveViewOrCreate('RenderView')
renderView.Background = [0.15, 0.15, 0.15]
Hide(reader, renderView)

# ── 1. Clip — expose vessel interior ─────────────────────────────────────────
clip = Clip(Input=reader)
clip.ClipType        = 'Plane'
clip.ClipType.Normal = [0.0, 0.0, 1.0]   # cut at z = 0
clip.ClipType.Origin = [0.0, 0.0, 0.0]
UpdatePipeline()

clipDisplay = Show(clip, renderView)
ColorBy(clipDisplay, ('CELLS', 'U', 'X'))
clipDisplay.RescaleTransferFunctionToDataRange(True)
clipDisplay.SetScalarBarVisibility(renderView, True)
clipDisplay.Opacity = 0.5

# ── 2. Slice — longitudinal midplane (z = 0) ─────────────────────────────────
slice1 = Slice(Input=reader)
slice1.SliceType        = 'Plane'
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.SliceType.Origin = [0.0, 0.0, 0.0]
UpdatePipeline()

sliceDisplay = Show(slice1, renderView)
ColorBy(sliceDisplay, ('CELLS', 'U', 'X'))
sliceDisplay.RescaleTransferFunctionToDataRange(True)

# ── 3. Stream Tracer — primary flow pattern visualization ─────────────────────
# Seed line placed near the inlet spanning 80% of the vessel diameter.
# In straight tubes: parallel lines.
# In junction/graft cases: closed loops reveal recirculation automatically.
streamTracer = StreamTracer(Input=reader, SeedType='Line')
streamTracer.Vectors                 = ['CELLS', 'U']
streamTracer.MaximumStreamlineLength = (bounds[1] - bounds[0]) * 2.0
streamTracer.SeedType.Point1         = [bounds[0] + 0.001, -RADIUS * 0.8, 0.0]
streamTracer.SeedType.Point2         = [bounds[0] + 0.001,  RADIUS * 0.8, 0.0]
streamTracer.SeedType.Resolution     = 20
UpdatePipeline()

streamDisplay = Show(streamTracer, renderView)
ColorBy(streamDisplay, ('POINTS', 'U', 'Magnitude'))
streamDisplay.RescaleTransferFunctionToDataRange(True)

# ── 4. Glyph — velocity arrows on the slice ───────────────────────────────────
glyph = Glyph(Input=slice1, GlyphType='Arrow')
glyph.OrientationArray = ['CELLS', 'U']
glyph.ScaleArray        = ['CELLS', 'U']
glyph.ScaleFactor       = RADIUS * 4.0    # auto-scales to vessel size
glyph.GlyphMode         = 'Every Nth Point'
glyph.Stride            = 2
UpdatePipeline()

glyphDisplay = Show(glyph, renderView)
ColorBy(glyphDisplay, ('POINTS', 'U', 'Magnitude'))

# ── 5. Plot Over Line — U_X profile chart near the outlet ────────────────────
plotLine = PlotOverLine(Input=reader)
plotLine.Point1     = [x_outlet, -RADIUS, 0.0]
plotLine.Point2     = [x_outlet,  RADIUS, 0.0]
plotLine.Resolution = 1000
UpdatePipeline()

chartView   = CreateView('XYChartView')
plotDisplay = Show(plotLine, chartView)
# Show only U_X — the axial velocity component (the Poiseuille parabola)
plotDisplay.SeriesVisibility = ['U_X']

# ── Reset camera and render ───────────────────────────────────────────────────
SetActiveView(renderView)
renderView.ResetCamera()
Render()

print(f"[standard_viz] Done — pipeline ready for: {CASE_NAME}\n")

# In pvpython the process would exit immediately without an event loop.
# Interact() blocks until the window is closed.
# Inside the ParaView GUI this call is skipped — the GUI has its own event loop.
if _IS_PVPYTHON:
    Interact()
