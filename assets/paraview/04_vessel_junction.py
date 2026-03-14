#!/usr/bin/env pvpython
"""
ParaView visualization -- Experiment 03: Womersley Pulsatile Flow Through a Vessel Junction
============================================================================================
Five-panel layout:
  - Left  60 % : 3-D render at peak systole (clip + slice + streamlines)
  - Right 40 % : 2×2 grid
      top-left     Outlet U_X profile  (x ≈ 95 mm, r2 = 3 mm, bars)
      top-right    Q(t) volumetric flow rate  [mL/s]
      bottom-left  Inlet U_X profile   (x ≈ 1 mm,  r1 = 5 mm, bars + fit)
      bottom-right Post-junction U_X   (x = 65 mm, r2 = 3 mm, bars) — shows recirculation

Geometry (Experiment 03):
  Upstream  : x = 0      → 0.06 m  radius r1 = 0.005 m  (donor vessel)
  Junction  : x = 0.06 m           sudden step contraction
  Downstream: x = 0.06   → 0.10 m  radius r2 = 0.003 m  (recipient vessel)

Q(t) and pressure are read from postProcessing files written by the
functionObjects block in system/controlDict (flowRateInlet, flowRateOutlet,
pressureProbes at x = 30 mm and x = 80 mm).

Pressure reconstruction (linear extrapolation within each uniform section):
  Probe 0 at x = 0.030 m (in r1 section)  ->  p0_kin [m2/s2]
  Probe 1 at x = 0.080 m (in r2 section)  ->  p1_kin [m2/s2]
  Each probe is close to its local outlet/inlet, so pressures are plotted
  directly (no HP extrapolation) to compare the two sections.
  Convert kinematic -> mmHg: p_mmHg = p_kin * rho / 133.322

How to run
----------
Option A -- pvpython:
    pvpython .../assets/paraview/04_vessel_junction.py \\
        --case .../run/04_vessel_junction
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
    "~/Rheology-Simulation-of-Vein-Grafts/run/04_vessel_junction"
)
R1       = 0.005    # upstream (donor) vessel radius [m]
R2       = 0.003    # downstream (recipient) vessel radius [m]
L_JCT    = 0.060    # junction x-position [m]
L_TOTAL  = 0.100    # total vessel length [m]
T_PERIOD = 0.857    # cardiac period [s]  (70 bpm)
PHI_PEAK = 0.175    # fraction of period at peak systole

RHO        = 1060.0
PA_TO_MMHG = 1.0 / 133.322
P_FACTOR   = RHO * PA_TO_MMHG   # p_kin [m2/s2] -> mmHg  (~= 7.9506)

M3_TO_ML = 1e6   # m3/s -> mL/s

# x-positions for line-probe slices
X_INLET   = 0.001          # just inside inlet (upstream, r1 section)
X_JDOWN   = L_JCT + 0.005  # 5 mm downstream of step (recirculation zone)
X_OUTLET  = L_TOTAL - 0.005 # 5 mm from outlet (downstream, r2 section)

# -- CLI args (pvpython only) -------------------------------------------------
if _IS_PVPYTHON:
    try:
        import argparse
        p = argparse.ArgumentParser()
        p.add_argument('--case', default=CASE_DIR)
        a, _ = p.parse_known_args()
        CASE_DIR = os.path.expanduser(a.case)
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
print("\n[03] Loading: " + FOAM_FILE)
reader = OpenFOAMReader(FileName=FOAM_FILE)
reader.MeshRegions = ['internalMesh']
reader.CellArrays  = ['U', 'p']
UpdatePipeline()

b     = reader.GetDataInformation().GetBounds()
x_min = b[0];  x_max = b[1]

# -- Peak-systole time --------------------------------------------------------
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
all_times = list(reader.TimestepValues)
if not all_times:
    raise RuntimeError("No time steps -- simulation incomplete?")

t_end         = all_times[-1]
t_cycle_start = (t_end // T_PERIOD) * T_PERIOD
t_peak        = min(all_times, key=lambda t: abs(t - (t_cycle_start + PHI_PEAK * T_PERIOD)))

print("[03]   t_peak = {:.3f} s  ({} time steps, t_end={:.3f} s)".format(
    t_peak, len(all_times), t_end))

scene.AnimationTime = t_peak
UpdatePipeline()

# =============================================================================
# POST-PROCESSING FILE HELPERS
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
    """Animated ProgrammableSource: draws a vertical green line at current time."""
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
    """Apply per-series colour/line styling for multiple series at once."""
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


def _profile_view(title, radius_m, y_min, y_max):
    """XYChartView for a U_X radial profile.  x-axis spans full diameter in mm."""
    v = CreateView('XYChartView')
    v.ChartTitle               = title
    v.BottomAxisTitle          = 'y  [mm]'
    v.LeftAxisTitle            = 'U_X  [m/s]'
    v.LeftAxisUseCustomRange   = 1
    v.LeftAxisRangeMinimum     = y_min
    v.LeftAxisRangeMaximum     = y_max
    v.BottomAxisUseCustomRange = 1
    v.BottomAxisRangeMinimum   = 0.0
    v.BottomAxisRangeMaximum   = radius_m * 2 * 1000
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
# Clip + slice through z=0 plane; streamlines seeded along x-axis
# to capture recirculation downstream of the step.
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

# Streamlines: seed at inlet (r1 range) AND inside downstream tube (r2 range)
# to capture both the main jet and any recirculation corner vortex
st = StreamTracer(Input=reader, SeedType='Line')
st.Vectors                 = ['CELLS', 'U']
st.MaximumStreamlineLength = L_TOTAL * 3
st.SeedType.Point1         = [x_min + 0.001, -R1 * 0.9, 0]
st.SeedType.Point2         = [x_min + 0.001,  R1 * 0.9, 0]
st.SeedType.Resolution     = 24
UpdatePipeline()
stDisp = Show(st, renderView)
ColorBy(stDisp, ('POINTS', 'U', 'Magnitude'))
stDisp.RescaleTransferFunctionToDataRange(True)

# Arrow glyphs on the midplane slice
gl = Glyph(Input=sl, GlyphType='Arrow')
gl.OrientationArray = ['CELLS', 'U']
gl.ScaleArray        = ['CELLS', 'U']
gl.ScaleFactor       = R2 * 4
gl.GlyphMode         = 'Every Nth Point'
gl.Stride            = 3
UpdatePipeline()
glDisp = Show(gl, renderView)
ColorBy(glDisp, ('POINTS', 'U', 'Magnitude'))

# =============================================================================
# CHART A -- Outlet U_X profile  (top-left, r2 range, bars)
# x = X_OUTLET ≈ 95 mm  (5 mm from outlet, inside r2 tube)
# Peak centre velocity expected ≈ (r1/r2)² × inlet peak ≈ 2.78 × inlet
# =============================================================================
x_out_mm = round(X_OUTLET * 1000, 1)

pol_out = PlotOverLine(Input=reader)
pol_out.Point1     = [X_OUTLET, -R2, 0]
pol_out.Point2     = [X_OUTLET,  R2, 0]
pol_out.Resolution = 500
UpdatePipeline()

calc_out = Calculator(Input=pol_out)
calc_out.AttributeType   = 'Point Data'
calc_out.ResultArrayName = 'arc_length_mm'
calc_out.Function        = 'arc_length * 1000'
UpdatePipeline()

outletChart = _profile_view(
    'Outlet  x = {} mm  (r2 = 3 mm)'.format(x_out_mm),
    radius_m=R2, y_min=-0.5, y_max=3.0
)
d_out = Show(calc_out, outletChart)
d_out.UseIndexForXAxis = 0
d_out.XArrayName       = 'arc_length_mm'
d_out.SeriesVisibility = ['U_X']
_style(d_out, 'U_X', 'U_X outlet  [m/s]', [0.80, 0.15, 0.15], 'bar')

# =============================================================================
# CHART B -- Q(t) volumetric flow rate  (top-right)
# Reads postProcessing/flowRateInlet|Outlet/.../surfaceFieldValue.dat
# =============================================================================
qChart = None

try:
    in_file  = _find_dat(os.path.join(POST, 'flowRateInlet'),  'surfaceFieldValue.dat')
    out_file = _find_dat(os.path.join(POST, 'flowRateOutlet'), 'surfaceFieldValue.dat')
    in_rows  = _read_dat(in_file)
    out_rows = _read_dat(out_file)

    times = [r[0] for r in in_rows]
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

    cur_q = _make_time_cursor(times[0], times[-1], q_lo, q_hi)
    d_cq  = Show(cur_q, qChart)
    d_cq.UseIndexForXAxis = 0
    d_cq.XArrayName       = 'time'
    _style_multi(d_cq, [
        {'name': 't_cursor', 'color': (0.18, 0.62, 0.17), 'thickness': 2, 'label': 'time'},
    ])

    print("[03]   Q(t) chart built ({} time steps)".format(len(times)))

except Exception as err:
    print("[03]   Q(t) chart SKIPPED: {}".format(err))
    print("       -> Run simulation with functionObjects in controlDict first.")
    qChart = _time_view('Q(t) -- run sim with functionObjects', 'Q  [mL/s]',
                        y_min=0.0, y_max=50.0)

# =============================================================================
# CHART C -- Inlet U_X profile  (bottom-left, r1 range, bars + fit)
# x = X_INLET ≈ 1 mm  (just inside inlet, upstream r1 tube)
# The Womersley profile differs from parabolic — fit overlay shows the gap.
# =============================================================================
x_in_mm = round(X_INLET * 1000, 1)

pol_in = PlotOverLine(Input=reader)
pol_in.Point1     = [X_INLET, -R1, 0]
pol_in.Point2     = [X_INLET,  R1, 0]
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
    'Inlet  x = {} mm  (r1 = 5 mm, fit updates on play)'.format(x_in_mm),
    radius_m=R1, y_min=0.0, y_max=1.2
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
_style(d_in_fit, 'U_X', 'Parabolic fit (live)', [0.122, 0.467, 0.706], 'line')

# =============================================================================
# CHART D -- Post-junction U_X profile  (bottom-right, r2 range, bars)
# x = X_JDOWN = 65 mm  (5 mm downstream of step at 60 mm)
# Key feature: annular recirculation ring during diastole (Womersley + step effect)
# Negative U_X near the outer wall (r ≈ r2) indicates reversed flow.
# =============================================================================
x_jd_mm = round(X_JDOWN * 1000, 1)

pol_jd = PlotOverLine(Input=reader)
pol_jd.Point1     = [X_JDOWN, -R2, 0]
pol_jd.Point2     = [X_JDOWN,  R2, 0]
pol_jd.Resolution = 500
UpdatePipeline()

calc_jd = Calculator(Input=pol_jd)
calc_jd.AttributeType   = 'Point Data'
calc_jd.ResultArrayName = 'arc_length_mm'
calc_jd.Function        = 'arc_length * 1000'
UpdatePipeline()

jdownChart = _profile_view(
    'Post-junction  x = {} mm  (r2 = 3 mm)'.format(x_jd_mm),
    radius_m=R2, y_min=-0.5, y_max=3.0
)
d_jd = Show(calc_jd, jdownChart)
d_jd.UseIndexForXAxis = 0
d_jd.XArrayName       = 'arc_length_mm'
d_jd.SeriesVisibility = ['U_X']
_style(d_jd, 'U_X', 'U_X post-junction  [m/s]', [0.58, 0.20, 0.56], 'bar')

# =============================================================================
# LAYOUT  --  3-D render left (60 %), 2×2 grid right (40 %)
#
#   Root(0) --SplitH 0.60-->  left(1) render    right(2)
#   right(2) --SplitV 0.50--> top(5)            bottom(6)
#   top(5)   --SplitH 0.50--> A outlet(11)      B Q(t)(12)
#   bottom(6) --SplitH 0.50-> C inlet(13)       D post-jct(14)
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
layout.AssignView(14, jdownChart)

# -- Final render -------------------------------------------------------------
SetActiveView(renderView)
renderView.ResetCamera()
Render()

print("\n[03] Done -- " + CASE_NAME)
print("  t_peak    = {:.3f} s  (peak systole of last cycle)".format(t_peak))
print("  Chart A   : Outlet U_X       (top-left)  -- red bars, r2=3mm, y in [-0.5, 3.0]")
print("  Chart B   : Q(t)             (top-right) -- inlet teal, outlet orange, green cursor")
print("  Chart C   : Inlet U_X        (bot-left)  -- red bars + live parabolic fit, r1=5mm")
print("  Chart D   : Post-junction UX (bot-right) -- purple bars, x=65mm, recirculation visible\n")

if _IS_PVPYTHON:
    Interact()
