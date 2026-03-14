#!/usr/bin/env pvpython
"""
ParaView visualization — Experiment 01: Simple Laminar Flow
============================================================
Pipeline for experiment 01_simple_laminar (steady Poiseuille flow in a
straight cylindrical vessel):
  1. Clip         — remove top half to expose the vessel interior
  2. Slice        — longitudinal midplane coloured by axial velocity (U_X)
  3. StreamTracer — streamlines coloured by U magnitude
  4. Glyph        — velocity arrows on the slice
  5. PlotOverLine — U_X spatial profile near the outlet (parabola validation)
  6. PlotOverLine — U_X spatial profile near the inlet  (1 mm from inlet face)

Layout (three panels):
  ┌─────────────────────┬──────────────────────┐
  │                     │  Outlet profile       │
  │   3-D Render        │  U_X vs y  [mm]       │
  │   (RenderView)      │  (PlotOverLine)       │
  │                     ├──────────────────────┤
  │                     │  Inlet profile        │
  │                     │  U_X vs y  [mm]       │
  │                     │  (PlotOverLine)       │
  └─────────────────────┴──────────────────────┘
          60 %                   40 %

How to run
----------
Option A — pvpython (command line, outside ParaView):
    pvpython ~/Rheology-Simulation-of-Vein-Grafts/assets/paraview/01_simple_laminar.py \\
        --case ~/Rheology-Simulation-of-Vein-Grafts/run/01_simple_laminar \\
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
print(f"\n[01_simple_laminar] Loading: {FOAM_FILE}")
reader = OpenFOAMReader(FileName=FOAM_FILE)
reader.MeshRegions = ['internalMesh']
reader.CellArrays  = ['U', 'p']
UpdatePipeline()

bounds = reader.GetDataInformation().GetBounds()
# bounds = (x_min, x_max, y_min, y_max, z_min, z_max)

# Inlet profile: 1 mm from the inlet face
x_inlet  = bounds[0] + 0.0
# Outlet profile: 90 % along the vessel toward the outlet
x_outlet = bounds[0] + 0.9 * (bounds[1] - bounds[0])

_vessel_len_mm = (bounds[1] - bounds[0]) * 1000
inlet_x_mm     = round((x_inlet  - bounds[0]) * 1000, 1)
outlet_x_mm    = round((x_outlet - bounds[0]) * 1000, 1)

print(f"[01_simple_laminar]   Vessel radius  : {RADIUS} m")
print(f"[01_simple_laminar]   Inlet profile  : x = {x_inlet:.4f} m  ({inlet_x_mm} mm)")
print(f"[01_simple_laminar]   Outlet profile : x = {x_outlet:.4f} m  ({outlet_x_mm} mm)")

# Jump to last time step (fully developed flow at t = 10 s)
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

# ── Helper: build a profile chart ────────────────────────────────────────────
def _make_profile_chart(reader, x_pos, y_min, y_max, title, fit_parabola=True):
    """Build an XYChartView for a U_X vs y [mm] cross-section profile.

    fit_parabola=True  (outlet):
        Fetches 200-point probe data, fits a degree-2 polynomial (physically
        correct for the Hagen-Poiseuille parabola), injects smooth 500-point
        curve via ProgrammableSource, overlays 25-point circle markers.

    fit_parabola=False (inlet / plug flow):
        Uses a 50-point PlotOverLine directly — data points connected by
        straight lines with circle markers, no curve fitting applied.
        The plug-flow profile is flat and non-parabolic; fitting a parabola
        would be physically wrong.
    """
    import numpy as np
    from paraview import servermanager as _sm

    BLUE = ['U_X', '0.122', '0.467', '0.706']
    LOG  = '[01_simple_laminar]'

    if fit_parabola:
        # ── High-res probe for polynomial fitting ─────────────────────────────
        pol_raw = PlotOverLine(Input=reader)
        pol_raw.Point1     = [x_pos, -RADIUS, 0.0]
        pol_raw.Point2     = [x_pos,  RADIUS, 0.0]
        pol_raw.Resolution = 200
        UpdatePipeline()

        # ── Low-res probe for markers (~1 per radial cell) ────────────────────
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

        # ── Degree-2 polynomial fit (parabola) ────────────────────────────────
        smooth_src = None
        try:
            raw_data = _sm.Fetch(pol_raw)
            n        = raw_data.GetNumberOfPoints()
            al_vtk   = raw_data.GetPointData().GetArray('arc_length')
            u_vtk    = raw_data.GetPointData().GetArray('U')
            y_arr    = np.array([al_vtk.GetValue(i) * 1000 for i in range(n)])
            ux_arr   = np.array([u_vtk.GetTuple3(i)[0]     for i in range(n)])

            coeffs  = np.polyfit(y_arr, ux_arr, 2)
            y_fine  = np.linspace(y_arr.min(), y_arr.max(), 500)
            ux_fine = np.polyval(coeffs, y_fine)

            smooth_src = ProgrammableSource()
            smooth_src.OutputDataSetType = 'vtkTable'
            smooth_src.Script = (
                "import vtk\n"
                f"y_data  = {y_fine.tolist()}\n"
                f"ux_data = {ux_fine.tolist()}\n"
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
            print(f"{LOG} Note: parabola fit failed ({ex}); showing raw curve")
            smooth_src = None

    else:
        # ── Direct line through data points (plug-flow / inlet) ───────────────
        # 50-point probe gives enough points to show the plug-flow shape clearly
        # without curve fitting. Data points are connected by straight lines.
        pol_pts = PlotOverLine(Input=reader)
        pol_pts.Point1     = [x_pos, -RADIUS, 0.0]
        pol_pts.Point2     = [x_pos,  RADIUS, 0.0]
        pol_pts.Resolution = 50
        UpdatePipeline()

        calc_pts = Calculator(Input=pol_pts)
        calc_pts.AttributeType   = 'Point Data'
        calc_pts.ResultArrayName = 'arc_length_mm'
        calc_pts.Function        = 'arc_length * 1000'
        UpdatePipeline()

        smooth_src = None   # no separate smooth source; calc_pts is the line

    # ── Chart view ────────────────────────────────────────────────────────────
    chart = CreateView('XYChartView')
    chart.ChartTitle               = title
    chart.BottomAxisTitle          = 'y  [mm]'
    chart.LeftAxisTitle            = 'U_X  [m/s]'
    chart.LeftAxisUseCustomRange   = 1
    chart.LeftAxisRangeMinimum     = y_min
    chart.LeftAxisRangeMaximum     = y_max
    chart.BottomAxisUseCustomRange = 1
    chart.BottomAxisRangeMinimum   = 0.0
    chart.BottomAxisRangeMaximum   = RADIUS * 2 * 1000

    # Layer 1: parabola line (fit_parabola) or raw connected line (plug flow)
    line_src  = smooth_src if smooth_src is not None else calc_pts
    disp_line = Show(line_src, chart)
    disp_line.UseIndexForXAxis = 0
    disp_line.XArrayName       = 'arc_length_mm'
    disp_line.SeriesVisibility = ['U_X']
    try:
        disp_line.SeriesColor         = BLUE
        disp_line.SeriesLineStyle     = ['U_X', '1']
        disp_line.SeriesLineThickness = ['U_X', '2']
        disp_line.SeriesMarkerStyle   = ['U_X', '0' if fit_parabola else '4']
        if not fit_parabola:
            disp_line.SeriesMarkerSize = ['U_X', '6']
    except Exception as e:
        print(f"{LOG} Note: line styling not applied ({e})")

    # Layer 2: separate circle markers (only when parabola line is shown)
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

# ── 5. PlotOverLine — outlet U_X profile (top-right panel) ───────────────────
# Parabola fit: H-P profile at fully-developed outlet is exactly parabolic.
outletChart = _make_profile_chart(
    reader, x_outlet, 0.0, 0.22,
    f'Outlet  ·  x = {outlet_x_mm} mm',
    fit_parabola=True
)

# ── 6. PlotOverLine — inlet U_X profile (bottom-right panel) ─────────────────
# No fitting: inlet is plug flow (flat top), parabola fit would be wrong.
# 50 data points connected by straight lines show the shape faithfully.
inletChart = _make_profile_chart(
    reader, x_inlet, 0.0, 0.22,
    f'Inlet  ·  x = {inlet_x_mm} mm',
    fit_parabola=False
)

# ── Layout ────────────────────────────────────────────────────────────────────
# Layout cell indices (ParaView binary-tree rule: children of cell i → 2i+1, 2i+2):
#   SplitHorizontal(0) → left=1, right=2
#   SplitVertical(2)   → top-right=5, bottom-right=6
layout = GetLayout(renderView)
layout.SplitHorizontal(0, 0.6)   # root → left (cell 1) | right (cell 2)
layout.SplitVertical(2, 0.5)     # right → top-right (cell 5) | bottom-right (cell 6)

layout.AssignView(1, renderView)   # 3-D render     (left)
layout.AssignView(5, outletChart)  # outlet profile (top-right)
layout.AssignView(6, inletChart)   # inlet  profile (bottom-right)

# ── Reset camera and render ───────────────────────────────────────────────────
SetActiveView(renderView)
renderView.ResetCamera()
Render()

print(f"[01_simple_laminar] Done — pipeline ready for: {CASE_NAME}\n")
print(f"  Top-right   : outlet profile (x = {outlet_x_mm} mm) — expect parabola peaking at ~0.19 m/s")
print(f"  Bottom-right: inlet  profile (x = {inlet_x_mm} mm)  — expect near-uniform plug flow\n")

if _IS_PVPYTHON:
    Interact()
