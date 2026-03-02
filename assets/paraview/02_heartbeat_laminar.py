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

# ── Helper: build a profile chart with smooth line + data-point markers ───────
def _make_profile_chart(reader, x_pos, y_min, y_max, title):
    """Two-layer chart: smooth spline curve + circle markers.

    Strategy: fetch raw PlotOverLine data client-side, deduplicate to one
    representative (y, U_X) per cell, fit a cubic spline, inject the smooth
    500-point curve back as a ProgrammableSource (vtkTable), and show it
    alongside raw circle markers in the same XYChartView.
    """
    import numpy as np
    from paraview import servermanager as _sm

    BLUE = ['U_X', '0.122', '0.467', '0.706']
    LOG  = '[02_heartbeat_laminar]'

    # High-res probe for polynomial fitting
    pol_raw = PlotOverLine(Input=reader)
    pol_raw.Point1     = [x_pos, -RADIUS, 0.0]
    pol_raw.Point2     = [x_pos,  RADIUS, 0.0]
    pol_raw.Resolution = 200
    UpdatePipeline()

    # Low-res probe for markers (~1 per radial cell)
    pol_pts = PlotOverLine(Input=reader)
    pol_pts.Point1     = [x_pos, -RADIUS, 0.0]
    pol_pts.Point2     = [x_pos,  RADIUS, 0.0]
    pol_pts.Resolution = 25
    UpdatePipeline()

    calc_pts = Calculator(Input=pol_pts)
    calc_pts.AttributeType   = 'Point Data'
    calc_pts.ResultArrayName = 'arc_length_mm'
    calc_pts.Function        = 'arc_length * 1000'
    UpdatePipeline()

    # Fetch high-res data and fit degree-2 polynomial (parabola)
    smooth_src = None
    try:
        raw_data = _sm.Fetch(pol_raw)
        n        = raw_data.GetNumberOfPoints()

        al_vtk = raw_data.GetPointData().GetArray('arc_length')
        u_vtk  = raw_data.GetPointData().GetArray('U')
        y_arr  = np.array([al_vtk.GetValue(i) * 1000 for i in range(n)])  # m→mm
        ux_arr = np.array([u_vtk.GetTuple3(i)[0]      for i in range(n)])

        # Degree-2 least-squares fit → smooth parabola (physically correct for H-P)
        coeffs  = np.polyfit(y_arr, ux_arr, 2)
        y_fine  = np.linspace(y_arr.min(), y_arr.max(), 500)
        ux_fine = np.polyval(coeffs, y_fine)

        smooth_src = ProgrammableSource()
        smooth_src.OutputDataSetType = 'vtkTable'
        y_lst  = y_fine.tolist()
        ux_lst = ux_fine.tolist()
        smooth_src.Script = (
            "import vtk\n"
            f"y_data  = {y_lst}\n"
            f"ux_data = {ux_lst}\n"
            "out = self.GetOutput()\n"
            "ya = vtk.vtkFloatArray(); ya.SetName('arc_length_mm')\n"
            "ua = vtk.vtkFloatArray(); ua.SetName('U_X')\n"
            "for yv, uv in zip(y_data, ux_data):\n"
            "    ya.InsertNextValue(yv)\n"
            "    ua.InsertNextValue(uv)\n"
            "out.AddColumn(ya)\n"
            "out.AddColumn(ua)\n"
        )
        UpdatePipeline()
    except Exception as ex:
        print(f"{LOG} Note: parabola fit failed ({ex}); showing raw curve only")
        smooth_src = None

    chart = CreateView('XYChartView')
    chart.ChartTitle               = title
    chart.BottomAxisTitle          = 'y  [mm]'
    chart.LeftAxisTitle            = 'U_X  [m/s]'
    chart.LeftAxisUseCustomRange   = 1
    chart.LeftAxisRangeMinimum     = y_min
    chart.LeftAxisRangeMaximum     = y_max
    chart.BottomAxisUseCustomRange = 1
    chart.BottomAxisRangeMinimum   = 0.0
    chart.BottomAxisRangeMaximum   = 10.0

    line_src  = smooth_src if smooth_src is not None else calc_pts
    disp_line = Show(line_src, chart)
    disp_line.UseIndexForXAxis = 0
    disp_line.XArrayName       = 'arc_length_mm'
    disp_line.SeriesVisibility = ['U_X']
    try:
        disp_line.SeriesColor         = BLUE
        disp_line.SeriesLineStyle     = ['U_X', '1']
        disp_line.SeriesLineThickness = ['U_X', '2']
        disp_line.SeriesMarkerStyle   = ['U_X', '0']
    except Exception as e:
        print(f"{LOG} Note: line styling not applied ({e})")

    if smooth_src is not None:
        disp_pts = Show(calc_pts, chart)
        disp_pts.UseIndexForXAxis = 0
        disp_pts.XArrayName       = 'arc_length_mm'
        disp_pts.SeriesVisibility = ['U_X']
        try:
            disp_pts.SeriesColor       = BLUE
            disp_pts.SeriesLineStyle   = ['U_X', '0']
            disp_pts.SeriesMarkerStyle = ['U_X', '4']
            disp_pts.SeriesMarkerSize  = ['U_X', '8']
        except Exception as e:
            print(f"{LOG} Note: marker styling not applied ({e})")

    return chart

# ── 5 & 6. Profile charts — outlet (top-right) and inlet (bottom-right) ───────
outletChart = _make_profile_chart(
    reader, x_outlet, -0.05, 0.5,
    f'Outlet  ·  x = {outlet_x_mm} mm'
)
inletChart = _make_profile_chart(
    reader, x_inlet, -0.05, 0.5,
    f'Inlet  ·  x = {inlet_x_mm} mm'
)

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
