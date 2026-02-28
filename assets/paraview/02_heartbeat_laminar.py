#!/usr/bin/env pvpython
"""
ParaView visualization — Experiment 02: Heartbeat Laminar Flow
==============================================================
Pipeline for experiment 02_heartbeat_laminar (pulsatile Poiseuille flow):
  1. Clip         — remove top half to expose the vessel interior
  2. StreamTracer — streamlines coloured by U magnitude
  3. Glyph        — velocity arrows on the longitudinal slice
  4. PlotOverLine — U_X spatial profile at the outlet (parabola check at
                    the current animation time step)
  5. ResampleWithDataset + PlotDataOverTime — U_X time history at the
                    outlet centreline (confirms heartbeat waveform and
                    periodic steady-state convergence)

Layout (three panels):
  ┌─────────────────────┬──────────────────────┐
  │                     │  Spatial profile      │
  │   3-D Render        │  U_X vs y  [m]        │
  │   (RenderView)      │  (PlotOverLine)       │
  │                     ├──────────────────────┤
  │                     │  Time history         │
  │                     │  U_X at centreline    │
  │                     │  (PlotDataOverTime)   │
  └─────────────────────┴──────────────────────┘
          60 %                   40 %

How to run
----------
Option A — pvpython (command line, outside ParaView):
    pvpython ~/Rheology-Simulation-of-Vein-Grafts/assets/paraview/02_heartbeat_laminar.py \\
        --case ~/Rheology-Simulation-of-Vein-Grafts/run/02_heartbeat_laminar \\
        --radius 0.005

Option B — ParaView Python Shell (Tools → Python Shell):
    Open the script and click Run Script.

Option C — ParaView Macro (one click):
    Tools → Macros → Add new macro → select this file → OK.

Parameters
----------
--case    Path to the run case folder. Must contain <case_name>.foam.
--radius  Inner vessel radius in metres (default: 0.005).
"""

import sys
import os

_IS_PVPYTHON = 'pvpython' in sys.executable

# ── Macro defaults ────────────────────────────────────────────────────────────
CASE_DIR = os.path.expanduser(
    "~/Rheology-Simulation-of-Vein-Grafts/run/02_heartbeat_laminar"
)
RADIUS = 0.005  # vessel inner radius [m]

# ── Parse command-line arguments (pvpython only) ──────────────────────────────
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

# ── Reset: close all open views and sources without saving ────────────────────
# Gives a clean slate every time the macro is run, regardless of what was
# previously open in the session.
try:
    ResetSession()
except Exception:
    try:
        for _proxy in list(GetSources().values()):
            Delete(_proxy)
        for _view in GetViews():
            Delete(_view)
    except Exception:
        pass

# ── Load the OpenFOAM case ────────────────────────────────────────────────────
print(f"\n[02_heartbeat_laminar] Loading: {FOAM_FILE}")
reader = OpenFOAMReader(FileName=FOAM_FILE)
reader.MeshRegions = ['internalMesh']
reader.CellArrays  = ['U', 'p']
UpdatePipeline()

bounds   = reader.GetDataInformation().GetBounds()
x_inlet  = bounds[0] + 0.01 * (bounds[1] - bounds[0])  # 1 % from inlet
x_outlet = bounds[0] + 0.9 * (bounds[1] - bounds[0])  # 90 % toward outlet

# Position labels used in chart titles
_vessel_len_mm = (bounds[1] - bounds[0]) * 1000
inlet_x_mm     = round(0.01 * _vessel_len_mm, 1)   # mm from vessel start
outlet_x_mm    = round(0.9 * _vessel_len_mm, 1)   # mm from vessel start

print(f"[02_heartbeat_laminar]   Vessel radius  : {RADIUS} m")
print(f"[02_heartbeat_laminar]   Inlet profile  : x = {x_inlet:.4f} m  ({inlet_x_mm} mm)")
print(f"[02_heartbeat_laminar]   Outlet profile : x = {x_outlet:.4f} m  ({outlet_x_mm} mm)")

# Animate to last time step (peak of last cardiac cycle)
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
clip.ClipType.Normal = [0.0, 0.0, 1.0]
clip.ClipType.Origin = [0.0, 0.0, 0.0]
UpdatePipeline()

clipDisplay = Show(clip, renderView)
ColorBy(clipDisplay, ('CELLS', 'U', 'X'))
clipDisplay.RescaleTransferFunctionToDataRange(True)
clipDisplay.SetScalarBarVisibility(renderView, True)
clipDisplay.Opacity = 0.5

# ── 2. Slice — longitudinal midplane for Glyph input ─────────────────────────
slice1 = Slice(Input=reader)
slice1.SliceType        = 'Plane'
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.SliceType.Origin = [0.0, 0.0, 0.0]
UpdatePipeline()

sliceDisplay = Show(slice1, renderView)
ColorBy(sliceDisplay, ('CELLS', 'U', 'X'))
sliceDisplay.RescaleTransferFunctionToDataRange(True)

# ── 3. StreamTracer — reveal flow pattern ────────────────────────────────────
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
glyph.ScaleFactor       = RADIUS * 4.0
glyph.GlyphMode         = 'Every Nth Point'
glyph.Stride            = 2
UpdatePipeline()

glyphDisplay = Show(glyph, renderView)
ColorBy(glyphDisplay, ('POINTS', 'U', 'Magnitude'))

# ── 5. PlotOverLine — U_X spatial profile at outlet (top-right panel) ─────────
# Shows the Hagen-Poiseuille parabola at the current time step.
# Animate the time slider to see how the profile amplitude pulses.
plotLine = PlotOverLine(Input=reader)
plotLine.Point1     = [x_outlet, -RADIUS, 0.0]
plotLine.Point2     = [x_outlet,  RADIUS, 0.0]
plotLine.Resolution = 1000
UpdatePipeline()

# Convert arc_length (m) → mm so X axis shows millimetres
calcMM = Calculator(Input=plotLine)
calcMM.AttributeType   = 'Point Data'
calcMM.ResultArrayName = 'arc_length_mm'
calcMM.Function        = 'arc_length * 1000'
UpdatePipeline()

def _make_profile_chart(calc_source, title):
    """Create a fixed-range XYChartView for a U_X vs y [mm] profile."""
    chart = CreateView('XYChartView')
    chart.ChartTitle          = title
    chart.BottomAxisTitle     = 'y  [mm]'
    chart.LeftAxisTitle       = 'U_X  [m/s]'
    chart.LeftAxisUseCustomRange   = 1
    chart.LeftAxisRangeMinimum     = -0.05
    chart.LeftAxisRangeMaximum     = 0.5
    chart.BottomAxisUseCustomRange = 1
    chart.BottomAxisRangeMinimum   = 0.0
    chart.BottomAxisRangeMaximum   = 10.0
    disp = Show(calc_source, chart)
    disp.UseIndexForXAxis = 0
    disp.XArrayName       = 'arc_length_mm'
    disp.SeriesVisibility = ['U_X']
    return chart

outletChart = _make_profile_chart(calcMM, f'Outlet  ·  x = {outlet_x_mm} mm')

# ── 6. PlotOverLine — U_X spatial profile at inlet (bottom-right panel) ───────
# Shows how the Hagen-Poiseuille parabola looks near the inlet at the current
# time step. Animating the time slider compares inlet and outlet amplitudes.
plotLineInlet = PlotOverLine(Input=reader)
plotLineInlet.Point1     = [x_inlet, -RADIUS, 0.0]
plotLineInlet.Point2     = [x_inlet,  RADIUS, 0.0]
plotLineInlet.Resolution = 1000
UpdatePipeline()

calcMM_inlet = Calculator(Input=plotLineInlet)
calcMM_inlet.AttributeType   = 'Point Data'
calcMM_inlet.ResultArrayName = 'arc_length_mm'
calcMM_inlet.Function        = 'arc_length * 1000'
UpdatePipeline()

inletChart = _make_profile_chart(calcMM_inlet, f'Inlet  ·  x = {inlet_x_mm} mm')

# ── Layout ────────────────────────────────────────────────────────────────────
# Layout cell indices (ParaView binary-tree rule: children of cell i → 2i+1, 2i+2):
#   SplitHorizontal(0) → left=1, right=2
#   SplitVertical(2)   → top-right=5, bottom-right=6
layout = GetLayout(renderView)
layout.SplitHorizontal(0, 0.6)   # root → left (cell 1) | right (cell 2)
layout.SplitVertical(2, 0.5)     # right → top-right (cell 5) | bottom-right (cell 6)

layout.AssignView(1, renderView)   # 3-D render (left)
layout.AssignView(5, outletChart)  # outlet profile (top-right)
layout.AssignView(6, inletChart)   # inlet  profile (bottom-right)

# ── Reset camera and render ───────────────────────────────────────────────────
SetActiveView(renderView)
renderView.ResetCamera()
Render()

print(f"[02_heartbeat_laminar] Done — pipeline ready for: {CASE_NAME}\n")
print("  Tip: Use the time slider (▶ Play) to animate the spatial profile")
print("       through the cardiac cycle and observe the pulsatile parabola.\n")

if _IS_PVPYTHON:
    Interact()
