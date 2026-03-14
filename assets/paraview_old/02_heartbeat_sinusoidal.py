#!/usr/bin/env pvpython
"""
ParaView visualization — Experiment 02: Heartbeat Laminar Flow
==============================================================
Pipeline for experiment 02_heartbeat_sinusoidal_no_initial_velocity (pulsatile Poiseuille flow):
  1. Clip         — remove top half to expose the vessel interior
  2. StreamTracer — streamlines coloured by U magnitude
  3. Glyph        — velocity arrows on the longitudinal slice
  4. PlotOverLine — U_X spatial profile at the outlet (animated parabola check)
  5. PlotOverLine — U_X spatial profile at the inlet  (animated parabola check)
  6. Flow-rate chart   — Q_inlet vs Q_outlet vs analytical Q(t) [mL/s]
  7. Inlet-pressure chart — p_inlet_sim vs p_outlet_sim [mmHg]

Layout (five panels):
  ┌──────────────┬──────────────────┬──────────────────┐
  │              │  Outlet profile  │  Flow rate Q(t)  │
  │  3-D Render  │  U_X vs y [mm]   │  [mL/s]          │
  │  (Render-    ├──────────────────┼──────────────────┤
  │   View)      │  Inlet  profile  │  Inlet/outlet p  │
  │              │  U_X vs y [mm]   │  p(t)  [mmHg]    │
  └──────────────┴──────────────────┴──────────────────┘
       35 %              32.5 %             32.5 %

Charts 6 & 7 read from postProcessing/ written by the controlDict functions
block (flowRateInlet, flowRateOutlet, pressureProbes). If those files are
absent the macro falls back to the original three-panel layout.

How to run
----------
Option A — pvpython (command line, outside ParaView):
    pvpython ~/Rheology-Simulation-of-Vein-Grafts/assets/paraview/02_heartbeat_sinusoidal_no_initial_velocity.py \\
        --case ~/Rheology-Simulation-of-Vein-Grafts/run/02_heartbeat_sinusoidal_no_initial_velocity \\
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
import math

_IS_PVPYTHON = 'pvpython' in sys.executable

# ── Macro defaults ────────────────────────────────────────────────────────────
CASE_DIR = os.path.expanduser(
    "~/Rheology-Simulation-of-Vein-Grafts/run/02_heartbeat_sinusoidal"
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

# ── Physics constants (must match 0/U and constant/transportProperties) ───────
NU        = 3.3e-6   # kinematic viscosity [m²/s]
RHO       = 1060.0   # blood density [kg/m³]
T_PERIOD  = 0.857    # cardiac period [s]  (70 bpm)
U_AVG     = 0.275    # mean centreline velocity [m/s]  (= (0.5+0.05)/2)
U_AMP     = 0.225    # waveform amplitude   [m/s]  (= (0.5-0.05)/2)
VESSEL_L  = 0.1      # tube length [m]

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
print(f"\n[02_heartbeat_sinusoidal_no_initial_velocity] Loading: {FOAM_FILE}")
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

print(f"[02_heartbeat_sinusoidal_no_initial_velocity]   Vessel radius  : {RADIUS} m")
print(f"[02_heartbeat_sinusoidal_no_initial_velocity]   Inlet profile  : x = {x_inlet:.4f} m  ({inlet_x_mm} mm)")
print(f"[02_heartbeat_sinusoidal_no_initial_velocity]   Outlet profile : x = {x_outlet:.4f} m  ({outlet_x_mm} mm)")

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
    LOG  = '[02_heartbeat_sinusoidal_no_initial_velocity]'

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

# ── Helpers for postProcessing time-series charts ─────────────────────────────

def _read_dat(filepath):
    """Parse a whitespace-delimited OpenFOAM postProcessing .dat file.

    Lines beginning with '#' or empty lines are skipped.
    Returns a list of rows, each row a list of floats.
    """
    rows = []
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            rows.append([float(x) for x in s.split()])
    return rows


def _make_table_source(columns):
    """Create a ProgrammableSource that outputs a vtkTable.

    columns — list of (name: str, values: list[float])
    """
    src = ProgrammableSource()
    src.OutputDataSetType = 'vtkTable'
    lines = ['import vtk', 'out = self.GetOutput()']
    for i, (name, vals) in enumerate(columns):
        lines.append(f'a{i} = vtk.vtkFloatArray(); a{i}.SetName({repr(name)})')
        lines.append(f'for v in {vals}: a{i}.InsertNextValue(v)')
        lines.append(f'out.AddColumn(a{i})')
    src.Script = '\n'.join(lines)
    UpdatePipeline()
    return src


def _style_series(disp, styles):
    """Apply per-series colour/line/marker styling to a chart display.

    styles — list of dicts with keys: name, color (R,G,B 0-1), line (0-4),
             thickness (int), marker (0=none 4=circle).
    """
    vis, color, lstyle, lthick, mstyle = [], [], [], [], []
    for s in styles:
        n = s['name']
        r, g, b = s.get('color', (0.2, 0.2, 0.2))
        vis   += [n]
        color += [n, str(r), str(g), str(b)]
        lstyle  += [n, str(s.get('line',      1))]
        lthick  += [n, str(s.get('thickness', 2))]
        mstyle  += [n, str(s.get('marker',    0))]
    try:
        disp.SeriesVisibility     = vis
        disp.SeriesColor          = color
        disp.SeriesLineStyle      = lstyle
        disp.SeriesLineThickness  = lthick
        disp.SeriesMarkerStyle    = mstyle
    except Exception as e:
        print(f'[02_heartbeat_sinusoidal_no_initial_velocity] Note: series styling skipped ({e})')


def _make_time_cursor(t_min, t_max, y_lo, y_hi, series_name='t_cursor'):
    """Time-aware ProgrammableSource that renders a green vertical line at the
    current animation time in an XYChartView.

    The source registers a TIME_RANGE so ParaView queries it at each animation
    step; each query returns exactly two rows (t_current, y_lo) and
    (t_current, y_hi), which the chart draws as a one-segment vertical line.
    """
    src = ProgrammableSource()
    src.OutputDataSetType = 'vtkTable'
    src.ScriptRequestInformation = (
        "import vtk\n"
        "ex = self.GetExecutive()\n"
        "oi = ex.GetOutputInformation(0)\n"
        "oi.Remove(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE())\n"
        f"oi.Append(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(), {t_min!r})\n"
        f"oi.Append(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(), {t_max!r})\n"
    )
    src.Script = (
        "import vtk\n"
        "ex = self.GetExecutive()\n"
        "oi = ex.GetOutputInformation(0)\n"
        "t = oi.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())\n"
        "out = self.GetOutput()\n"
        f"xa = vtk.vtkFloatArray(); xa.SetName('time')\n"
        f"ya = vtk.vtkFloatArray(); ya.SetName({series_name!r})\n"
        f"xa.InsertNextValue(t); ya.InsertNextValue({y_lo!r})\n"
        f"xa.InsertNextValue(t); ya.InsertNextValue({y_hi!r})\n"
        "out.AddColumn(xa)\n"
        "out.AddColumn(ya)\n"
    )
    UpdatePipeline()
    return src


# ── 7. Flow-rate conservation chart ──────────────────────────────────────────
# Reads postProcessing/flowRateInlet(Outlet)/<time>/surfaceFieldValue.dat.
# Plots |Q_inlet| and Q_outlet alongside the analytical sinusoidal waveform.
# Units converted to mL/s (1 m³/s = 1e6 mL/s) for readability.

import glob as _glob
import traceback as _tb

LOG = '[02_heartbeat_sinusoidal_no_initial_velocity]'
POST = os.path.join(CASE_DIR, 'postProcessing')
M3_TO_ML = 1e6   # m³/s → mL/s

print(f'{LOG} postProcessing folder: {POST}')
print(f'{LOG} Exists: {os.path.isdir(POST)}')

flowRateChart  = None
pressureChart  = None


def _find_dat(base_dir, filename):
    """Return the first matching file found under base_dir/<any_time>/<filename>.

    OpenFOAM may use '0', '0.04', or the startTime value as the time-directory
    name depending on version and writeControl settings.
    """
    pattern = os.path.join(base_dir, '*', filename)
    matches = sorted(_glob.glob(pattern))
    if matches:
        print(f'{LOG}   Found: {matches[0]}')
        return matches[0]
    raise FileNotFoundError(f'No file matching {pattern}')


try:
    inlet_file  = _find_dat(os.path.join(POST, 'flowRateInlet'),  'surfaceFieldValue.dat')
    outlet_file = _find_dat(os.path.join(POST, 'flowRateOutlet'), 'surfaceFieldValue.dat')
    inlet_rows  = _read_dat(inlet_file)
    outlet_rows = _read_dat(outlet_file)

    times = [r[0] for r in inlet_rows]
    q_in  = [-r[1] * M3_TO_ML for r in inlet_rows]   # negate: inlet phi < 0
    q_out = [ r[1] * M3_TO_ML for r in outlet_rows]

    q_lo = min(min(q_in), min(q_out)) * 0.95
    q_hi = max(max(q_in), max(q_out)) * 1.05

    tbl_q = _make_table_source([
        ('time',     times),
        ('Q_inlet',  q_in),
        ('Q_outlet', q_out),
    ])

    flowRateChart = CreateView('XYChartView')
    flowRateChart.ChartTitle      = 'Flow Rate Conservation  Q(t)'
    flowRateChart.BottomAxisTitle = 'Time  [s]'
    flowRateChart.LeftAxisTitle   = 'Q  [mL/s]'

    disp_q = Show(tbl_q, flowRateChart)
    disp_q.UseIndexForXAxis = 0
    disp_q.XArrayName       = 'time'
    _style_series(disp_q, [
        {'name': 'Q_inlet',  'color': (0.122, 0.467, 0.706), 'line': 1, 'thickness': 2},
        {'name': 'Q_outlet', 'color': (1.0,   0.498, 0.055), 'line': 1, 'thickness': 2},
    ])

    _cursor_q = _make_time_cursor(times[0], times[-1], q_lo, q_hi)
    _disp_cq  = Show(_cursor_q, flowRateChart)
    _disp_cq.UseIndexForXAxis = 0
    _disp_cq.XArrayName       = 'time'
    _style_series(_disp_cq, [
        {'name': 't_cursor', 'color': (0.18, 0.62, 0.17), 'line': 1, 'thickness': 2, 'marker': 0},
    ])
    print(f'{LOG} Flow-rate chart built ({len(times)} time steps)')

except Exception:
    print(f'{LOG} Flow-rate chart FAILED:')
    _tb.print_exc()


# ── 8. Inlet / outlet pressure chart ─────────────────────────────────────────
# Reads postProcessing/pressureProbes/<time>/p which contains kinematic pressure
# at two interior probe locations:
#   Probe 0  (x = 0.02375 m, 10th axial cell ≈ quarter-length)
#   Probe 1  (x = 0.07375 m, 30th axial cell ≈ three-quarter-length)
#
# HP linear profile: p(x) = p_inlet · (1 − x/L)
#   ⟹  p_inlet  = 1.475 · p₀ − 0.475 · p₁   (extrapolate to x = 0)
#       p_outlet = −0.525 · p₀ + 1.525 · p₁  (extrapolate to x = L)
#
# Units: kinematic p → Pa  (×ρ),  Pa → mmHg  (÷133.322).

PA_TO_MMHG = 1.0 / 133.322   # 1 Pa = 0.007501 mmHg

try:
    _probe_file = _find_dat(os.path.join(POST, 'pressureProbes'), 'p')
    _probe_rows = _read_dat(_probe_file)

    _ts_p   = [r[0] for r in _probe_rows]
    _p0_raw = [r[1] for r in _probe_rows]   # probe 0  (x = 0.02375 m)
    _p1_raw = [r[2] for r in _probe_rows]   # probe 1  (x = 0.07375 m)

    # HP extrapolation to x = 0 (inlet) and x = L (outlet), converted to mmHg
    _p_inlet_mmhg  = [(1.475 * p0 - 0.475 * p1) * RHO * PA_TO_MMHG
                      for p0, p1 in zip(_p0_raw, _p1_raw)]
    _p_outlet_mmhg = [(-0.525 * p0 + 1.525 * p1) * RHO * PA_TO_MMHG
                      for p0, p1 in zip(_p0_raw, _p1_raw)]

    _p_hi = max(max(_p_inlet_mmhg), max(_p_outlet_mmhg)) * 1.05

    _tbl_p = _make_table_source([
        ('time',         _ts_p),
        ('p_inlet_sim',  _p_inlet_mmhg),
        ('p_outlet_sim', _p_outlet_mmhg),
    ])

    pressureChart = CreateView('XYChartView')
    pressureChart.ChartTitle             = 'Inlet / Outlet Pressure  p(t)'
    pressureChart.BottomAxisTitle        = 'Time  [s]'
    pressureChart.LeftAxisTitle          = 'p  [mmHg]'
    pressureChart.LeftAxisUseCustomRange = 1
    pressureChart.LeftAxisRangeMinimum   = 0.0
    pressureChart.LeftAxisRangeMaximum   = _p_hi

    _disp_p = Show(_tbl_p, pressureChart)
    _disp_p.UseIndexForXAxis = 0
    _disp_p.XArrayName       = 'time'
    _style_series(_disp_p, [
        {'name': 'p_inlet_sim',  'color': (0.122, 0.467, 0.706), 'line': 1, 'thickness': 2},
        {'name': 'p_outlet_sim', 'color': (0.839, 0.153, 0.157), 'line': 1, 'thickness': 2},
    ])

    _cursor_p = _make_time_cursor(_ts_p[0], _ts_p[-1], 0.0, _p_hi)
    _disp_cp  = Show(_cursor_p, pressureChart)
    _disp_cp.UseIndexForXAxis = 0
    _disp_cp.XArrayName       = 'time'
    _style_series(_disp_cp, [
        {'name': 't_cursor', 'color': (0.18, 0.62, 0.17), 'line': 1, 'thickness': 2, 'marker': 0},
    ])
    print(f'{LOG} Inlet/outlet-pressure chart built ({len(_ts_p)} time steps)')

except Exception:
    print(f'{LOG} Pressure chart FAILED:')
    _tb.print_exc()


# ── Layout ────────────────────────────────────────────────────────────────────
# Binary-tree cell indices: children of cell i → 2i+1 (left/top), 2i+2 (right/bottom)
#
# Five-panel (both extra charts):
#   Split H(0,0.35)→1,2  H(2,0.5)→5,6  V(5,0.5)→11,12  V(6,0.5)→13,14
#   1=render  11=outlet  12=inlet  13=flowRate  14=pressure
#
# Four-panel (one extra chart — flow rate only):
#   Split H(0,0.35)→1,2  V(2,0.5)→5,6  H(5,0.5)→11,12
#   1=render  11=outlet  12=flowRate  6=inlet
#
# Three-panel fallback:
#   Split H(0,0.6)→1,2  V(2,0.5)→5,6
#   1=render  5=outlet  6=inlet

layout = GetLayout(renderView)
print(f'{LOG} flowRateChart={flowRateChart is not None}  pressureChart={pressureChart is not None}')

if flowRateChart and pressureChart:
    # ── Five-panel ──
    layout.SplitHorizontal(0, 0.35)
    layout.SplitHorizontal(2, 0.5)
    layout.SplitVertical(5, 0.5)
    layout.SplitVertical(6, 0.5)
    layout.AssignView(1,  renderView)
    layout.AssignView(11, outletChart)
    layout.AssignView(12, inletChart)
    layout.AssignView(13, flowRateChart)
    layout.AssignView(14, pressureChart)
    print(f'{LOG} Five-panel layout assigned.')

elif flowRateChart:
    # ── Four-panel (flow rate only) ──
    layout.SplitHorizontal(0, 0.35)
    layout.SplitVertical(2, 0.5)
    layout.SplitHorizontal(5, 0.5)
    layout.AssignView(1,  renderView)
    layout.AssignView(11, outletChart)
    layout.AssignView(12, flowRateChart)
    layout.AssignView(6,  inletChart)
    print(f'{LOG} Four-panel layout assigned (pressure chart missing).')

else:
    # ── Three-panel fallback ──
    layout.SplitHorizontal(0, 0.6)
    layout.SplitVertical(2, 0.5)
    layout.AssignView(1, renderView)
    layout.AssignView(5, outletChart)
    layout.AssignView(6, inletChart)
    print(f'{LOG} Three-panel fallback layout (no postProcessing charts available).')


# ── Reset camera and render ───────────────────────────────────────────────────
SetActiveView(renderView)
renderView.ResetCamera()
Render()

print(f'\n{LOG} Done — pipeline ready for: {CASE_NAME}')
print('  Tip: Use the time slider (▶ Play) to animate the spatial profile')
print('       through the cardiac cycle and observe the pulsatile parabola.')
print('  Flow-rate chart  : blue=Q_inlet  orange=Q_outlet  green=time cursor')
print('  Pressure chart   : blue=p_inlet_sim  red=p_outlet_sim  green=time cursor  [mmHg]')

if _IS_PVPYTHON:
    Interact()
