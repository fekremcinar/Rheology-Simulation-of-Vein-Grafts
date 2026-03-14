#!/usr/bin/env pvpython
"""
ParaView visualization -- Experiment 02: Heartbeat Laminar Flow
===============================================================
Five-panel layout:
  - Left 60 %  : 3-D render at peak systole
  - Right 40 % : 2x2 grid
      top-left     Outlet U_X profile  (bars only)
      top-right    Q(t) volumetric flow rate  [mL/s]
      bottom-left  Inlet U_X profile   (bars + live parabola fit)
      bottom-right Inlet/outlet pressure p(t) [mmHg]

Q(t) and p(t) are read from OpenFOAM postProcessing files written by the
functionObjects block in system/controlDict (flowRateInlet, flowRateOutlet,
pressureProbes).  The simulation MUST be run with those functionObjects
enabled -- see controlDict.  If the postProcessing folder is absent the
script still runs but those two charts will be skipped.

Pressure reconstruction (Hagen-Poiseuille linear extrapolation):
  Probe 0 at x = 0.025 m (x/L = 0.25)  ->  p0_kin [m2/s2]
  Probe 1 at x = 0.075 m (x/L = 0.75)  ->  p1_kin [m2/s2]
  p_inlet  =  1.5 * p0 - 0.5 * p1  (extrapolate to x = 0)
  p_outlet = -0.5 * p0 + 1.5 * p1  (extrapolate to x = L)
  Convert kinematic -> mmHg: p_mmHg = p_kin * rho / 133.322  (~= p_kin * 7.9506)

How to run
----------
Option A -- pvpython:
    pvpython .../assets/paraview/02_biphasic_heartbeat.py \\
        --case .../run/02_biphasic_heartbeat --radius 0.005
Option B -- ParaView Python Shell: open & Run Script.
Option C -- ParaView Macro: Tools -> Macros -> Add new macro.
"""

import sys
import os
import math
import glob as _glob

_IS_PVPYTHON = 'pvpython' in sys.executable

# -- Macro defaults -----------------------------------------------------------
CASE_DIR = os.path.expanduser(
    "~/Rheology-Simulation-of-Vein-Grafts/run/02_biphasic_heartbeat"
)
RADIUS   = 0.0005   # vessel inner radius [m]
T_PERIOD = 0.857    # cardiac period [s]  (70 bpm)
PHI_PEAK = 0.175    # fraction of period at peak systole

RHO        = 1060.0
PA_TO_MMHG = 1.0 / 133.322
P_FACTOR   = RHO * PA_TO_MMHG   # p_kin [m2/s2] -> mmHg  (~= 7.9506)

# HP extrapolation coefficients for probes at x=0.025, x=0.075, L=0.1
# p_inlet  = C0_IN  * p0 + C1_IN  * p1
# p_outlet = C0_OUT * p0 + C1_OUT * p1
C0_IN  =  1.5;  C1_IN  = -0.5
C0_OUT = -0.5;  C1_OUT =  1.5

M3_TO_ML = 1e6   # m3/s -> mL/s

# -- CLI args (pvpython only) -------------------------------------------------
if _IS_PVPYTHON:
    try:
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument('--case',   default=CASE_DIR)
        p.add_argument('--radius', type=float, default=RADIUS)
        a, _ = p.parse_known_args()
        CASE_DIR = os.path.expanduser(a.case)
        RADIUS   = a.radius
    except Exception:
        pass

CASE_NAME = os.path.basename(CASE_DIR.rstrip('/'))
FOAM_FILE = os.path.join(CASE_DIR, CASE_NAME + '.foam')
POST      = os.path.join(CASE_DIR, 'postProcessing')

# -- ParaView -----------------------------------------------------------------
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

try:
    ResetSession()
except Exception:
    try:
        for _px in list(GetSources().values()):
            Delete(_px)
        for _vw in GetViews():
            Delete(_vw)
    except Exception:
        pass

# -- Load case ----------------------------------------------------------------
print("\n[02] Loading: " + FOAM_FILE)
reader = OpenFOAMReader(FileName=FOAM_FILE)
reader.MeshRegions = ['internalMesh']
reader.CellArrays  = ['U', 'p']
UpdatePipeline()

b     = reader.GetDataInformation().GetBounds()
x_min = b[0];  x_max = b[1];  L_vessel = x_max - x_min

x_inlet  = x_min + 0.001
x_outlet = x_min + 0.9 * L_vessel

# -- Peak-systole time --------------------------------------------------------
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
all_times = list(reader.TimestepValues)
if not all_times:
    raise RuntimeError("No time steps -- simulation incomplete?")

t_end         = all_times[-1]
t_cycle_start = (t_end // T_PERIOD) * T_PERIOD
t_peak        = min(all_times, key=lambda t: abs(t - (t_cycle_start + PHI_PEAK * T_PERIOD)))

print("[02]   t_peak = {:.3f} s  ({} time steps, t_end={:.3f} s)".format(
    t_peak, len(all_times), t_end))

scene.AnimationTime = t_peak
UpdatePipeline()

# =============================================================================
# POST-PROCESSING FILE HELPERS
# Ported from assets/paraview_old/02_heartbeat_sinusoidal.py
# =============================================================================
def _find_dat(base_dir, filename):
    """Return path of first matching postProcessing file (any time subdir)."""
    pattern = os.path.join(base_dir, '*', filename)
    matches = sorted(_glob.glob(pattern))
    if matches:
        return matches[0]
    raise FileNotFoundError('No file matching ' + pattern)


def _read_dat(filepath):
    """Parse OpenFOAM postProcessing .dat file into list of float rows."""
    rows = []
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            rows.append([float(x) for x in s.split()])
    return rows


def _make_table_source(columns):
    """ProgrammableSource -> vtkTable with named float columns.

    columns: list of (name: str, values: list[float])
    """
    src = ProgrammableSource()
    src.OutputDataSetType = 'vtkTable'
    lines = ['import vtk', 'out = self.GetOutput()']
    for i, (name, vals) in enumerate(columns):
        lines.append('a{} = vtk.vtkFloatArray(); a{}.SetName({!r})'.format(i, i, name))
        lines.append('for v in {}: a{}.InsertNextValue(v)'.format(vals, i))
        lines.append('out.AddColumn(a{})'.format(i))
    src.Script = '\n'.join(lines)
    UpdatePipeline()
    return src


def _make_time_cursor(t_min, t_max, y_lo, y_hi, series_name='t_cursor'):
    """Animated ProgrammableSource: draws a vertical green line at current time.

    Uses ScriptRequestInformation to register a TIME_RANGE so ParaView
    queries this source at every animation step and redraws the cursor.
    """
    src = ProgrammableSource()
    src.OutputDataSetType = 'vtkTable'
    src.ScriptRequestInformation = (
        "import vtk\n"
        "ex = self.GetExecutive()\n"
        "oi = ex.GetOutputInformation(0)\n"
        "oi.Remove(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE())\n"
        "oi.Append(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(), {!r})\n"
        "oi.Append(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(), {!r})\n"
    ).format(t_min, t_max)
    src.Script = (
        "import vtk\n"
        "ex = self.GetExecutive()\n"
        "oi = ex.GetOutputInformation(0)\n"
        "t = oi.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())\n"
        "out = self.GetOutput()\n"
        "xa = vtk.vtkFloatArray(); xa.SetName('time')\n"
        "ya = vtk.vtkFloatArray(); ya.SetName({!r})\n"
        "xa.InsertNextValue(t); ya.InsertNextValue({!r})\n"
        "xa.InsertNextValue(t); ya.InsertNextValue({!r})\n"
        "out.AddColumn(xa)\n"
        "out.AddColumn(ya)\n"
    ).format(series_name, float(y_lo), float(y_hi))
    UpdatePipeline()
    return src


# =============================================================================
# STYLING HELPERS
# =============================================================================
def _style(disp, arr, label, rgb, chart_type='line'):
    """Label, colour, chart-type for a single series."""
    try:
        ct = '0' if chart_type == 'line' else '2'
        disp.SeriesChartType     = [arr, ct]
        disp.SeriesColor         = [arr, str(rgb[0]), str(rgb[1]), str(rgb[2])]
        disp.SeriesLineStyle     = [arr, '1' if chart_type == 'line' else '0']
        disp.SeriesLineThickness = [arr, '2']
        disp.SeriesMarkerStyle   = [arr, '0']
        disp.SeriesLabel         = [arr, label]
    except Exception as ex:
        print("    _style ({}): {}".format(arr, ex))


def _style_multi(disp, styles):
    """Apply per-series colour/line styling for multiple series at once.

    styles: list of dicts with keys: name, color (R,G,B), line, thickness
    """
    vis, color, lstyle, lthick, mstyle, labels = [], [], [], [], [], []
    for s in styles:
        n   = s['name']
        r, g, b = s.get('color', (0.2, 0.2, 0.2))
        vis    += [n]
        color  += [n, str(r), str(g), str(b)]
        lstyle += [n, str(s.get('line',      1))]
        lthick += [n, str(s.get('thickness', 2))]
        mstyle += [n, str(s.get('marker',    0))]
        labels += [n, s.get('label', n)]
    try:
        disp.SeriesVisibility    = vis
        disp.SeriesColor         = color
        disp.SeriesLineStyle     = lstyle
        disp.SeriesLineThickness = lthick
        disp.SeriesMarkerStyle   = mstyle
        disp.SeriesLabel         = labels
    except Exception as ex:
        print("    _style_multi: {}".format(ex))


def _profile_view(title, y_min, y_max):
    """XYChartView for a U_X radial profile."""
    v = CreateView('XYChartView')
    v.ChartTitle               = title
    v.BottomAxisTitle          = 'y  [mm]'
    v.LeftAxisTitle            = 'U_X  [m/s]'
    v.LeftAxisUseCustomRange   = 1
    v.LeftAxisRangeMinimum     = y_min
    v.LeftAxisRangeMaximum     = y_max
    v.BottomAxisUseCustomRange = 1
    v.BottomAxisRangeMinimum   = 0.0
    v.BottomAxisRangeMaximum   = RADIUS * 2 * 1000
    return v


def _time_view(title, y_label, y_min, y_max):
    """XYChartView for a time-series f(t)."""
    v = CreateView('XYChartView')
    v.ChartTitle               = title
    v.BottomAxisTitle          = 't  [s]'
    v.LeftAxisTitle            = y_label
    v.BottomAxisUseCustomRange = 1
    v.BottomAxisRangeMinimum   = 0.0
    v.BottomAxisRangeMaximum   = t_end
    v.LeftAxisUseCustomRange   = 1
    v.LeftAxisRangeMinimum     = y_min
    v.LeftAxisRangeMaximum     = y_max
    return v


# =============================================================================
# 3-D RENDER VIEW
# =============================================================================
renderView = GetActiveViewOrCreate('RenderView')
renderView.Background = [0.15, 0.15, 0.15]
Hide(reader, renderView)

clip = Clip(Input=reader)
clip.ClipType.Normal = [0, 0, 1]
clip.ClipType.Origin = [0, 0, 0]
UpdatePipeline()
clipDisp = Show(clip, renderView)
ColorBy(clipDisp, ('CELLS', 'U', 'X'))
clipDisp.RescaleTransferFunctionToDataRange(True)
clipDisp.SetScalarBarVisibility(renderView, True)
clipDisp.Opacity = 0.5

sl = Slice(Input=reader)
sl.SliceType.Normal = [0, 0, 1]
sl.SliceType.Origin = [0, 0, 0]
UpdatePipeline()
slDisp = Show(sl, renderView)
ColorBy(slDisp, ('CELLS', 'U', 'X'))
slDisp.RescaleTransferFunctionToDataRange(True)

st = StreamTracer(Input=reader, SeedType='Line')
st.Vectors                 = ['CELLS', 'U']
st.MaximumStreamlineLength = L_vessel * 2
st.SeedType.Point1         = [x_min + 0.001, -RADIUS * 0.8, 0]
st.SeedType.Point2         = [x_min + 0.001,  RADIUS * 0.8, 0]
st.SeedType.Resolution     = 20
UpdatePipeline()
stDisp = Show(st, renderView)
ColorBy(stDisp, ('POINTS', 'U', 'Magnitude'))
stDisp.RescaleTransferFunctionToDataRange(True)

gl = Glyph(Input=sl, GlyphType='Arrow')
gl.OrientationArray = ['CELLS', 'U']
gl.ScaleArray        = ['CELLS', 'U']
gl.ScaleFactor       = RADIUS * 4
gl.GlyphMode         = 'Every Nth Point'
gl.Stride            = 2
UpdatePipeline()
glDisp = Show(gl, renderView)
ColorBy(glDisp, ('POINTS', 'U', 'Magnitude'))

# =============================================================================
# CHART A -- Outlet U_X profile  (top-left of 2x2, bars only)
# =============================================================================
x_out_mm = round((x_outlet - x_min) * 1000, 1)

pol_out = PlotOverLine(Input=reader)
pol_out.Point1     = [x_outlet, -RADIUS, 0]
pol_out.Point2     = [x_outlet,  RADIUS, 0]
pol_out.Resolution = 500
UpdatePipeline()

calc_out = Calculator(Input=pol_out)
calc_out.AttributeType   = 'Point Data'
calc_out.ResultArrayName = 'arc_length_mm'
calc_out.Function        = 'arc_length * 1000'
UpdatePipeline()

outletChart = _profile_view(
    'Outlet  x = {} mm'.format(x_out_mm),
    y_min=-0.2, y_max=1.0
)
d_out = Show(calc_out, outletChart)
d_out.UseIndexForXAxis = 0
d_out.XArrayName       = 'arc_length_mm'
d_out.SeriesVisibility = ['U_X']
_style(d_out, 'U_X', 'U_X outlet  [m/s]', [0.80, 0.15, 0.15], 'bar')

# =============================================================================
# CHART B -- Q(t) volumetric flow rate  (top-right of 2x2)
# Reads postProcessing/flowRateInlet|Outlet/.../surfaceFieldValue.dat
# =============================================================================
qChart        = None
flowRateChart = None

try:
    in_file  = _find_dat(os.path.join(POST, 'flowRateInlet'),  'surfaceFieldValue.dat')
    out_file = _find_dat(os.path.join(POST, 'flowRateOutlet'), 'surfaceFieldValue.dat')
    in_rows  = _read_dat(in_file)
    out_rows = _read_dat(out_file)

    times = [r[0] for r in in_rows]
    # phi at inlet is negative (outward face normal points -x); negate to get Q > 0
    q_in  = [-r[1] * M3_TO_ML for r in in_rows]
    q_out = [ r[1] * M3_TO_ML for r in out_rows]

    q_hi = max(max(q_in), max(q_out)) * 1.10
    q_lo = min(0.0, min(min(q_in), min(q_out)) * 0.95)

    tbl_q = _make_table_source([
        ('time',     times),
        ('Q_inlet',  q_in),
        ('Q_outlet', q_out),
    ])

    qChart = _time_view('Flow Rate  Q(t)', 'Q  [mL/s]', y_min=q_lo, y_max=q_hi)

    d_q = Show(tbl_q, qChart)
    d_q.UseIndexForXAxis = 0
    d_q.XArrayName       = 'time'
    _style_multi(d_q, [
        {'name': 'Q_inlet',  'color': (0.122, 0.467, 0.706), 'label': 'Q inlet  [mL/s]'},
        {'name': 'Q_outlet', 'color': (0.890, 0.467, 0.122), 'label': 'Q outlet [mL/s]'},
    ])

    # Animated vertical time-cursor (green)
    cur_q = _make_time_cursor(times[0], times[-1], q_lo, q_hi)
    d_cq  = Show(cur_q, qChart)
    d_cq.UseIndexForXAxis = 0
    d_cq.XArrayName       = 'time'
    _style_multi(d_cq, [
        {'name': 't_cursor', 'color': (0.18, 0.62, 0.17), 'thickness': 2, 'label': 'time'},
    ])

    print("[02]   Q(t) chart built ({} time steps)".format(len(times)))

except Exception as err:
    print("[02]   Q(t) chart SKIPPED: {}".format(err))
    print("       -> Run simulation with functionObjects in controlDict first.")
    qChart = _time_view('Q(t) -- run sim with functionObjects', 'Q  [mL/s]',
                        y_min=0.0, y_max=50.0)

# =============================================================================
# CHART C -- Inlet U_X profile  (bottom-left of 2x2, bars + live parabola fit)
# =============================================================================
x_in_mm = round((x_inlet - x_min) * 1000, 1)

pol_in = PlotOverLine(Input=reader)
pol_in.Point1     = [x_inlet, -RADIUS, 0]
pol_in.Point2     = [x_inlet,  RADIUS, 0]
pol_in.Resolution = 500
UpdatePipeline()

calc_in = Calculator(Input=pol_in)
calc_in.AttributeType   = 'Point Data'
calc_in.ResultArrayName = 'arc_length_mm'
calc_in.Function        = 'arc_length * 1000'
UpdatePipeline()

fit_in = ProgrammableFilter(Input=pol_in)
fit_in.OutputDataSetType = 'vtkTable'
fit_in.Script = (
    "import vtk, numpy as np\n"
    "inp = self.GetInputDataObject(0, 0)\n"
    "out = self.GetOutput()\n"
    "ya = vtk.vtkFloatArray(); ya.SetName('arc_length_mm')\n"
    "ua = vtk.vtkFloatArray(); ua.SetName('U_X')\n"
    "n = inp.GetNumberOfPoints()\n"
    "if n >= 3:\n"
    "    al = inp.GetPointData().GetArray('arc_length')\n"
    "    uv = inp.GetPointData().GetArray('U')\n"
    "    if al and uv:\n"
    "        y_arr  = np.array([al.GetValue(i)*1000 for i in range(n)])\n"
    "        ux_arr = np.array([uv.GetTuple3(i)[0]  for i in range(n)])\n"
    "        try:\n"
    "            c = np.polyfit(y_arr, ux_arr, 2)\n"
    "            y_fit  = np.linspace(y_arr.min(), y_arr.max(), 500)\n"
    "            ux_fit = np.polyval(c, y_fit)\n"
    "        except Exception:\n"
    "            y_fit, ux_fit = y_arr, ux_arr\n"
    "        for yv, fv in zip(y_fit, ux_fit):\n"
    "            ya.InsertNextValue(float(yv))\n"
    "            ua.InsertNextValue(float(fv))\n"
    "out.AddColumn(ya)\n"
    "out.AddColumn(ua)\n"
)
UpdatePipeline()

inletChart = _profile_view(
    'Inlet  x = {} mm  (fit updates on play)'.format(x_in_mm),
    y_min=0.0, y_max=1.2
)
d_in_bar = Show(calc_in, inletChart)
d_in_bar.UseIndexForXAxis = 0
d_in_bar.XArrayName       = 'arc_length_mm'
d_in_bar.SeriesVisibility = ['U_X']
_style(d_in_bar, 'U_X', 'Sampled (500 pts)', [0.80, 0.15, 0.15], 'bar')

d_in_fit = Show(fit_in, inletChart)
d_in_fit.UseIndexForXAxis = 0
d_in_fit.XArrayName       = 'arc_length_mm'
d_in_fit.SeriesVisibility = ['U_X']
_style(d_in_fit, 'U_X', 'Fitted parabola (live)', [0.122, 0.467, 0.706], 'line')

# =============================================================================
# CHART D -- Inlet/outlet pressure p(t) [mmHg]  (bottom-right of 2x2)
# Reads postProcessing/pressureProbes/.../p
# HP extrapolation: p_inlet = 1.5*p0 - 0.5*p1,  p_outlet = -0.5*p0 + 1.5*p1
# =============================================================================
pressChart = None

try:
    p_file = _find_dat(os.path.join(POST, 'pressureProbes'), 'p')
    p_rows = _read_dat(p_file)

    ts_p   = [r[0] for r in p_rows]
    p0_raw = [r[1] for r in p_rows]   # probe 0 at x = 0.025 m
    p1_raw = [r[2] for r in p_rows]   # probe 1 at x = 0.075 m

    # HP linear extrapolation, then kinematic -> mmHg
    p_in_mmhg  = [(C0_IN  * p0 + C1_IN  * p1) * P_FACTOR for p0, p1 in zip(p0_raw, p1_raw)]
    p_out_mmhg = [(C0_OUT * p0 + C1_OUT * p1) * P_FACTOR for p0, p1 in zip(p0_raw, p1_raw)]

    p_hi = max(max(p_in_mmhg), max(p_out_mmhg)) * 1.10
    p_hi = max(p_hi, 0.1)   # at least 0.1 mmHg on axis

    tbl_p = _make_table_source([
        ('time',         ts_p),
        ('p_inlet_mmHg',  p_in_mmhg),
        ('p_outlet_mmHg', p_out_mmhg),
    ])

    pressChart = _time_view('Pressure  p(t)', 'p  [mmHg]', y_min=0.0, y_max=p_hi)

    d_p = Show(tbl_p, pressChart)
    d_p.UseIndexForXAxis = 0
    d_p.XArrayName       = 'time'
    _style_multi(d_p, [
        {'name': 'p_inlet_mmHg',  'color': (0.122, 0.467, 0.706), 'label': 'p inlet  [mmHg]'},
        {'name': 'p_outlet_mmHg', 'color': (0.839, 0.153, 0.157), 'label': 'p outlet [mmHg]'},
    ])

    # Animated vertical time-cursor (green)
    cur_p = _make_time_cursor(ts_p[0], ts_p[-1], 0.0, p_hi)
    d_cp  = Show(cur_p, pressChart)
    d_cp.UseIndexForXAxis = 0
    d_cp.XArrayName       = 'time'
    _style_multi(d_cp, [
        {'name': 't_cursor', 'color': (0.18, 0.62, 0.17), 'thickness': 2, 'label': 'time'},
    ])

    print("[02]   p(t) chart built ({} time steps)".format(len(ts_p)))

except Exception as err:
    print("[02]   p(t) chart SKIPPED: {}".format(err))
    print("       -> Run simulation with functionObjects in controlDict first.")
    pressChart = _time_view('p(t) -- run sim with functionObjects', 'p  [mmHg]',
                            y_min=0.0, y_max=2.0)

# =============================================================================
# LAYOUT  --  2x2 grid on the right
#
#   Root(0) --SplitH 0.60-->  left(1) render    right(2)
#   right(2) --SplitV 0.50--> top(5)            bottom(6)
#   top(5)   --SplitH 0.50--> A outlet(11)      B Q(t)(12)
#   bottom(6) --SplitH 0.50-> C inlet(13)       D p(t)(14)
# =============================================================================
layout = GetLayout(renderView)
layout.SplitHorizontal(0, 0.6)
layout.SplitVertical(2, 0.5)
layout.SplitHorizontal(5, 0.5)
layout.SplitHorizontal(6, 0.5)

layout.AssignView(1,  renderView)
layout.AssignView(11, outletChart)
layout.AssignView(12, qChart)
layout.AssignView(13, inletChart)
layout.AssignView(14, pressChart)

# -- Final render -------------------------------------------------------------
SetActiveView(renderView)
renderView.ResetCamera()
Render()

print("\n[02] Done -- " + CASE_NAME)
print("  t_peak    = {:.3f} s  (peak systole of last cycle)".format(t_peak))
print("  Chart A   : Outlet U_X  (top-left)  -- red bars, y in [-0.2, 1.0]")
print("  Chart B   : Q(t)        (top-right) -- inlet teal, outlet orange, green cursor")
print("  Chart C   : Inlet U_X   (bot-left)  -- red bars + live parabola fit")
print("  Chart D   : p(t) mmHg   (bot-right) -- inlet blue, outlet red, green cursor\n")

if _IS_PVPYTHON:
    Interact()
