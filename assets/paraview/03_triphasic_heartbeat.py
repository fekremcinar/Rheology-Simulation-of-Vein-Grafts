#!/usr/bin/env pvpython
"""
ParaView visualization -- Experiment 03: Triphasic Heartbeat Laminar Flow
==========================================================================
Five-panel layout:
  - Left 60 %  : 3-D render at peak systole
  - Right 40 % : 2x2 grid
      top-left     Outlet U_X profile  (bars only)
      top-right    Q(t) volumetric flow rate  [mL/s]
      bottom-left  Inlet U_X profile   (bars + live parabola fit)
      bottom-right Hemodynamic metrics (dP, WSS, Re, RRT normalised by H-P)

Q(t) and hemodynamic data are read from OpenFOAM postProcessing files written
by the functionObjects block in system/controlDict (flowRateInlet, flowRateOutlet,
pAvgInlet, pAvgOutlet).  The simulation MUST be run with those functionObjects
enabled -- see controlDict.  If the postProcessing folder is absent the
script still runs but those charts will be skipped.

How to run
----------
Option A -- pvpython:
    pvpython .../assets/paraview/03_triphasic_heartbeat.py \\
        --case .../run/03_triphasic_heartbeat --radius 0.0005
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
    "~/Rheology-Simulation-of-Vein-Grafts/run/03_triphasic_heartbeat"
)
RADIUS   = 0.0005   # vessel inner radius [m]
T_PERIOD = 0.857    # cardiac period [s]  (70 bpm)
PHI_PEAK = 0.175    # fraction of period at peak systole

RHO        = 1060.0
MU         = 0.0035          # dynamic viscosity [Pa·s]
NU         = MU / RHO        # kinematic viscosity [m²/s]

A_INLET    = math.pi * RADIUS**2

# Mean flow parameters (triphasic waveform, Murray scaling to 1 mm artery)
Q_MEAN     = 1.656e-8        # [m³/s]
U_MEAN     = Q_MEAN / A_INLET

# Triphasic waveform peak/reversal extremes (Q_norm_peak=12, Q_norm_min=-3.2)
# U_center = 2 * Q_norm * U_MEAN  (parabolic profile)
U_CENTER_PEAK = 2.0 * 12.0  * U_MEAN   # ≈ +0.506 m/s (forward peak)
U_CENTER_MIN  = 2.0 * (-3.2) * U_MEAN  # ≈ -0.135 m/s (reversal)

# Hagen-Poiseuille reference values (Poiseuille at Q_mean)
L_TUBE     = 0.01            # vessel length [m]
DP_HP_PA   = 8.0 * MU * Q_MEAN * L_TUBE / (math.pi * RADIUS**4)   # [Pa]
WSS_HP     = 4.0 * MU * Q_MEAN / (math.pi * RADIUS**3)            # [Pa]
RE_THEORY  = U_MEAN * (2.0 * RADIUS) / NU
RRT_HP     = 1.0 / WSS_HP                                          # [Pa^-1]

M3_TO_ML   = 1e6   # m³/s -> mL/s

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

# NOTE: Q(t) can go negative during the flow reversal phase (~t=0.24-0.39 s).
# This is physically correct for the triphasic peripheral/femoral waveform.
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

x_inlet  = x_min
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


def _profile_view(title, y_min, y_max, x_max=1.0):
    """XYChartView for a U_X radial profile."""
    v = CreateView('XYChartView')
    v.ChartTitle               = title
    v.BottomAxisTitle          = 'y  [mm]'
    v.LeftAxisTitle            = 'U_X  [cm/s]'
    v.LeftAxisUseCustomRange   = 1
    v.LeftAxisRangeMinimum     = y_min
    v.LeftAxisRangeMaximum     = y_max
    v.BottomAxisUseCustomRange = 1
    v.BottomAxisRangeMinimum   = 0.0
    v.BottomAxisRangeMaximum   = x_max
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
clipDisp.SetScalarBarVisibility(renderView, True)
clipDisp.Opacity = 0.5

sl = Slice(Input=reader)
sl.SliceType.Normal = [0, 0, 1]
sl.SliceType.Origin = [0, 0, 0]
UpdatePipeline()
slDisp = Show(sl, renderView)
ColorBy(slDisp, ('CELLS', 'U', 'X'))

st = StreamTracer(Input=reader, SeedType='Line')
st.Vectors                 = ['CELLS', 'U']
st.MaximumStreamlineLength = L_vessel * 2
st.SeedType.Point1         = [x_min + 0.001, -RADIUS * 0.8, 0]
st.SeedType.Point2         = [x_min + 0.001,  RADIUS * 0.8, 0]
st.SeedType.Resolution     = 20
UpdatePipeline()
stDisp = Show(st, renderView)
ColorBy(stDisp, ('POINTS', 'U', 'Magnitude'))

gl = Glyph(Input=sl, GlyphType='Arrow')
gl.OrientationArray = ['CELLS', 'U']
gl.ScaleArray        = ['CELLS', 'U']
gl.ScaleFactor       = RADIUS * 4
gl.GlyphMode         = 'Every Nth Point'
gl.Stride            = 2
UpdatePipeline()
glDisp = Show(gl, renderView)
ColorBy(glDisp, ('POINTS', 'U', 'Magnitude'))

# Fix colormap range to cover the full pulsatile velocity range so the
# render doesn't saturate at the low-velocity state (RescaleToDataRange
# only captures the current timestep which may not be t_peak).
# Triphasic waveform: includes flow reversal → min is negative.
lut = GetColorTransferFunction('U')
lut.RescaleTransferFunction(U_CENTER_MIN, U_CENTER_PEAK)
pwf = GetOpacityTransferFunction('U')
pwf.RescaleTransferFunction(U_CENTER_MIN, U_CENTER_PEAK)

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

calc_out_cms = Calculator(Input=calc_out)
calc_out_cms.AttributeType   = 'Point Data'
calc_out_cms.ResultArrayName = 'U_X_cms'
calc_out_cms.Function        = 'U_X * 100'
UpdatePipeline()

outletChart = _profile_view(
    'Outlet  x = {} mm'.format(x_out_mm),
    y_min=-16, y_max=54, x_max=1.0
)
d_out = Show(calc_out_cms, outletChart)
d_out.UseIndexForXAxis = 0
d_out.XArrayName       = 'arc_length_mm'
d_out.SeriesVisibility = ['U_X_cms']
_style(d_out, 'U_X_cms', 'U_X outlet  [cm/s]', [0.80, 0.15, 0.15], 'bar')

# Parabola fit on outlet — values in cm/s
fit_out = ProgrammableFilter(Input=pol_out)
fit_out.OutputDataSetType = 'vtkTable'
fit_out.Script = (
    "import vtk, numpy as np\n"
    "inp = self.GetInputDataObject(0, 0)\n"
    "out = self.GetOutput()\n"
    "ya = vtk.vtkFloatArray(); ya.SetName('arc_length_mm')\n"
    "ua = vtk.vtkFloatArray(); ua.SetName('U_X_cms')\n"
    "n = inp.GetNumberOfPoints()\n"
    "if n >= 3:\n"
    "    al = inp.GetPointData().GetArray('arc_length')\n"
    "    uv = inp.GetPointData().GetArray('U')\n"
    "    if al and uv:\n"
    "        y_arr  = np.array([al.GetValue(i)*1000 for i in range(n)])\n"
    "        ux_arr = np.array([uv.GetTuple3(i)[0]*100 for i in range(n)])\n"
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
d_out_fit = Show(fit_out, outletChart)
d_out_fit.UseIndexForXAxis = 0
d_out_fit.XArrayName       = 'arc_length_mm'
d_out_fit.SeriesVisibility = ['U_X_cms']
_style(d_out_fit, 'U_X_cms', 'Fitted parabola (live)', [0.122, 0.467, 0.706], 'line')

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

    q_hi = max(max(q_in), max(q_out)) * 1.15
    q_lo = min(0.0, min(min(q_in), min(q_out))) * 1.15

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

calc_in_cms = Calculator(Input=calc_in)
calc_in_cms.AttributeType   = 'Point Data'
calc_in_cms.ResultArrayName = 'U_X_cms'
calc_in_cms.Function        = 'U_X * 100'
UpdatePipeline()

fit_in = ProgrammableFilter(Input=pol_in)
fit_in.OutputDataSetType = 'vtkTable'
fit_in.Script = (
    "import vtk, numpy as np\n"
    "inp = self.GetInputDataObject(0, 0)\n"
    "out = self.GetOutput()\n"
    "ya = vtk.vtkFloatArray(); ya.SetName('arc_length_mm')\n"
    "ua = vtk.vtkFloatArray(); ua.SetName('U_X_cms')\n"
    "n = inp.GetNumberOfPoints()\n"
    "if n >= 3:\n"
    "    al = inp.GetPointData().GetArray('arc_length')\n"
    "    uv = inp.GetPointData().GetArray('U')\n"
    "    if al and uv:\n"
    "        y_arr  = np.array([al.GetValue(i)*1000 for i in range(n)])\n"
    "        ux_arr = np.array([uv.GetTuple3(i)[0]*100  for i in range(n)])\n"
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
    'Inlet  x = {} mm'.format(x_in_mm),
    y_min=-16, y_max=54, x_max=1.0
)
d_in_bar = Show(calc_in_cms, inletChart)
d_in_bar.UseIndexForXAxis = 0
d_in_bar.XArrayName       = 'arc_length_mm'
d_in_bar.SeriesVisibility = ['U_X_cms']
_style(d_in_bar, 'U_X_cms', 'Sampled (500 pts)', [0.80, 0.15, 0.15], 'bar')

d_in_fit = Show(fit_in, inletChart)
d_in_fit.UseIndexForXAxis = 0
d_in_fit.XArrayName       = 'arc_length_mm'
d_in_fit.SeriesVisibility = ['U_X_cms']
_style(d_in_fit, 'U_X_cms', 'Fitted parabola (live)', [0.122, 0.467, 0.706], 'line')

# =============================================================================
# CHART D -- Hemodynamic Metrics  (bottom-right)
#
# Reads pAvgInlet/pAvgOutlet and flowRateInlet from postProcessing.
# Computes per time step:
#   dP(t)     = (pAvgInlet - pAvgOutlet) * rho          [Pa]
#   WSS_HP(t) = 4 * mu * Q(t) / (pi * r^3)              [Pa]
#   Re(t)     = Q(t) / A / (2r) / nu                    [-]
#   RRT(t)    = 1 / WSS_HP(t)                           [Pa^-1]  (OSI=0 steady)
#
# All metrics are normalised by their Hagen-Poiseuille reference value.
# A healthy fully-developed Poiseuille flow gives ratio = 1.0 for all metrics.
#
# Thrombosis thresholds (not normalised):
#   TAWSS < 0.4 Pa  -> thrombogenic wall region
#   RRT   > 20 Pa^-1 -> high platelet retention risk
# =============================================================================
hemoChart = None
try:
    in_file = _find_dat(os.path.join(POST, 'flowRateInlet'), 'surfaceFieldValue.dat')
    pi_file = _find_dat(os.path.join(POST, 'pAvgInlet'),     'surfaceFieldValue.dat')
    po_file = _find_dat(os.path.join(POST, 'pAvgOutlet'),    'surfaceFieldValue.dat')

    q_rows  = _read_dat(in_file)
    pi_rows = _read_dat(pi_file)
    po_rows = _read_dat(po_file)

    times_h = [r[0] for r in q_rows]

    # Q(t): inlet phi < 0 during forward flow (inward face normal) -> negate.
    # During triphasic reversal phi > 0 (outflow) -> q_vals < 0.
    q_vals = [-r[1] for r in q_rows]   # [m^3/s]

    # dP(t): kinematic -> physical [Pa]
    dp_vals = [(pi[1] - po[1]) * RHO for pi, po in zip(pi_rows, po_rows)]

    # WSS_HP(t): from Q using Hagen-Poiseuille relation tau_w = 4*mu*Q/(pi*r^3)
    wss_vals = [4.0 * MU * q / (math.pi * RADIUS**3) for q in q_vals]

    # Re(t) = U_mean * D / nu = Q / A * 2r / nu
    re_vals = [q / A_INLET * (2.0 * RADIUS) / NU for q in q_vals]

    # RRT(t) = 1 / WSS.
    # Guard: during triphasic flow reversal Q crosses zero → WSS → 0 → RRT → ∞.
    # Use a threshold of 1% of HP reference WSS to avoid division-by-near-zero spikes.
    rrt_vals = [1.0 / w if abs(w) > WSS_HP * 0.01 else 0.0 for w in wss_vals]

    # Normalise each series by its H-P reference value
    def _safe_div(a, b):
        return a / b if abs(b) > 1e-15 else 0.0

    dp_norm  = [_safe_div(v, DP_HP_PA)  for v in dp_vals]
    wss_norm = [_safe_div(v, WSS_HP)    for v in wss_vals]
    re_norm  = [_safe_div(v, RE_THEORY) for v in re_vals]
    rrt_norm = [_safe_div(v, RRT_HP)    for v in rrt_vals]
    ones     = [1.0] * len(times_h)

    # Clip all normalised values to avoid residual spikes at zero-crossings
    # dominating the y-axis scale. Physiological range for a healthy tube is
    # [-2, 15] (peak normalised Re ≈ Q_norm_peak = 12 for triphasic waveform).
    CLIP = 15.0
    dp_norm  = [max(-CLIP, min(CLIP, v)) for v in dp_norm]
    wss_norm = [max(-CLIP, min(CLIP, v)) for v in wss_norm]
    re_norm  = [max(-CLIP, min(CLIP, v)) for v in re_norm]
    rrt_norm = [max(-CLIP, min(CLIP, v)) for v in rrt_norm]

    tbl_h = _make_table_source([
        ('time',         times_h),
        ('dP_norm',      dp_norm),
        ('WSS_norm',     wss_norm),
        ('Re_norm',      re_norm),
        ('RRT_norm',     rrt_norm),
        ('HP_reference', ones),
    ])

    all_norm = dp_norm + wss_norm + re_norm + rrt_norm
    # Use only the last 30 % of time steps to set the y-axis range so the
    # startup transient does not dominate (data is already clipped).
    tail_start = max(1, int(len(times_h) * 0.70))
    tail_vals  = (dp_norm[tail_start:] + wss_norm[tail_start:] +
                  re_norm[tail_start:] + rrt_norm[tail_start:])
    if not tail_vals:
        tail_vals = all_norm if all_norm else [1.0]
    y_lo = min(0.0, min(tail_vals) * 1.1)
    y_hi = max(max(tail_vals) * 1.3, 1.5)

    hemoChart = _time_view('Hemodynamic Metrics', 'metric / HP value  [-]', y_min=y_lo, y_max=y_hi)

    d_h = Show(tbl_h, hemoChart)
    d_h.UseIndexForXAxis = 0
    d_h.XArrayName       = 'time'
    _style_multi(d_h, [
        {'name': 'dP_norm',
         'color': (0.122, 0.467, 0.706), 'thickness': 2,
         'label': 'dP / dP_HP  ({:.1f} Pa)'.format(DP_HP_PA)},
        {'name': 'WSS_norm',
         'color': (0.839, 0.153, 0.157), 'thickness': 2,
         'label': 'TAWSS / TAWSS_HP  ({:.3f} Pa)'.format(WSS_HP)},
        {'name': 'Re_norm',
         'color': (0.549, 0.337, 0.294), 'thickness': 2,
         'label': 'Re / Re_HP  ({:.0f})'.format(RE_THEORY)},
        {'name': 'RRT_norm',
         'color': (0.498, 0.498, 0.498), 'thickness': 2,
         'label': 'RRT / RRT_HP  ({:.1f} Pa^-1)'.format(RRT_HP)},
        {'name': 'HP_reference',
         'color': (0.18, 0.62, 0.17), 'thickness': 1, 'line': 2,
         'label': 'H-P reference (= 1.0)'},
    ])

    cur_h = _make_time_cursor(times_h[0], times_h[-1], y_lo, y_hi,
                              series_name='t_cursor_h')
    d_ch  = Show(cur_h, hemoChart)
    d_ch.UseIndexForXAxis = 0
    d_ch.XArrayName       = 'time'
    _style_multi(d_ch, [
        {'name': 't_cursor_h', 'color': (0.18, 0.62, 0.17), 'thickness': 2, 'label': 'current time'},
    ])

    print("[03]   Hemodynamic metrics chart built ({} time steps)".format(len(times_h)))
    print("[03]   H-P refs: dP={:.2f} Pa, WSS={:.4f} Pa, Re={:.0f}, RRT={:.2f} Pa^-1".format(
        DP_HP_PA, WSS_HP, RE_THEORY, RRT_HP))

except Exception as err:
    print("[03]   Hemodynamic metrics chart SKIPPED: {}".format(err))
    print("       -> Run simulation with functionObjects in controlDict first.")
    hemoChart = _time_view(
        'Hemodynamic Metrics -- run sim with functionObjects',
        'metric / HP value  [-]', y_min=0.0, y_max=2.0
    )

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
layout.AssignView(14, hemoChart)

# -- Final render -------------------------------------------------------------
SetActiveView(renderView)
renderView.ResetCamera()
Render()

print("\n[03] Done -- " + CASE_NAME)
print("  t_peak    = {:.3f} s  (peak systole of last cycle)".format(t_peak))
print("  Chart A   : Outlet U_X  (top-left)  -- red bars, y in [-16, 54] cm/s")
print("  Chart B   : Q(t)        (top-right) -- inlet teal, outlet orange, green cursor")
print("  Chart C   : Inlet U_X   (bot-left)  -- red bars + live parabola fit")
print("  Chart D   : Hemodynamics (bot-right) -- dP, WSS, Re, RRT normalised by H-P\n")

if _IS_PVPYTHON:
    Interact()
