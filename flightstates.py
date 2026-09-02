"""
flightstates.py — flight track in, state list out.

The unit is the SECOND, not the thermal. Every second of the track gets two
properties — is the pilot circling, and what is the height doing — and a
segment is a contiguous run of equal seconds. A thermal is therefore a run of
"circling and climbing", no more and no less, and it stands beside "straight
and climbing" as an equal.

Two kinds of segment, and only two — everything else is a number:

    K  circling      the pilot turns; this is where the wind is measured
    G  straight      everything else; carries w, the movement of the air

There are no classes. Every straight segment states w in m/s — the vertical
movement of the air, the glider's own sink at the flown airspeed already
removed. Whether the pilot climbed is the height difference of the segment;
how strongly the air rose is w. Classes would only round that number.

Line format (one flight, one line, fields separated by semicolons):

    id;glider;yyyymmdd;SEG|SEG|...;END

  SEG straight (G,T,L):   A,t,h,lat,lon
  SEG circling (K,k):     A,t,h,lat,lon,turns,drift_kmh,drift_from_deg
  END:                      t,h,lat,lon

  A       kind (one letter)
  t       second of the day UTC at which the segment BEGINS
  h       altitude MSL in m at the beginning
  lat,lon position at the beginning, five decimals (about 1 m, as in kk7)

The END of a segment is the BEGINNING of the next — which is why it is not
written twice. Only the last point of the flight needs an entry of its own.

The wind belongs to the circling and is written only there. A long climb
is cut into pieces of about 4 full turns (K_TURNS_PIECE), each with its own
drift — the wind profile through the climb, not one mean. Measured by the
method of thermal_strategy.py: at least 1.5 full turns, at least 30 s, four
fixes trimmed at each end (entry and exit are not circling), drift above
40 km/h discarded. If the conditions are not met, both wind fields stay empty;
the turn count is written anyway.

Dependencies: numpy, pandas, scipy. Nothing else — this file runs on its own.

    python3 flightstates.py flights/*.IGC > states.txt
    python3 flightstates.py --delta flights/*.IGC > states.txt   (compact)

The second field is the glider type, taken from the HFGTY record of the IGC
header. It names no person, and it is what a sink polar has to be fitted per:
estimated per flight that quantity scatters by 0.35 m/s — as much as the air
movement one wants to measure — while per glider type it settles to about
0.03 m/s. The line stores the type and nothing derived from it; the polar, the
airspeed and the correction all belong in the analysis, not in the file.

The first field is YOURS. This script never puts a name there; by default it
writes the single letter "P". Fill it with whatever suits you: an anonymised
pilot key, a flight number of your own, a constant, or nothing at all.

    --id TEXT     write TEXT in the first field (default "P")
    --id name     use the file name — only sensible on your own archive,
                  since IGC file names usually carry the pilot's name

If you put an anonymised PILOT key there, one thing becomes possible that no
single flight can show: how a pilot's technique develops over the years — when
the circles get tighter, when the straight climbing starts to appear, how the
glides change. Two things could come of that, and both are yours to decide. I
would study it, but ANONYMOUSLY ONLY — never a named pilot, in any publication
of mine. And you could turn it into something for your own users, a personal
progress report for paying XCTrack subscribers, for instance. Nothing here
needs the field to be filled.

Author: Jörg Korner, 2026. Free to use.
"""
import glob as _glob
import sys, re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

# ---------------------------------------------------------------- thresholds
TURN_WINDOW = 20      # s, centred window for the heading change
TURN_LIMIT  = 60.0    # deg per 20 s -> 3 deg/s sustained = circling
WIND_FACTOR  = 0.85    # ground speed must exceed 0.85 x wind, else it is drift
CLIMBS      = 0.10    # m/s, above this the piece counts as climbing
# The glider's own sink in calm air, MEASURED: the mode of the vario over
# 4.19 million seconds of straight flight (489 alpine cross-country flights,
# 1164 hours). The distribution has a single broad peak and no shoulders, so
# the class boundaries below are conventions — but they are anchored on that
# measurement and expressed in the movement of the AIR,
#     w = vario - CALM_SINK,
# which is the physical quantity, instead of the vario, which also depends on
# how fast the pilot chose to fly.
CALM_SINK   = -1.11   # m/s
W_LIFT      = 0.60    # m/s, above this the air is carrying   -> T
W_SINK      = -0.60   # m/s, below this the air is sinking    -> L
W_STRONG    = -1.50   # m/s, below this it is sinking hard    -> S
CARRIES     = CALM_SINK + W_LIFT   # -0.51 m/s, kept as the old name
MIN_RUN_S     = 20      # s, shorter circling/straight pieces are absorbed
# Simplification of the straight pieces (Douglas-Peucker). A corner is set
# only where a straight line would otherwise be off by more than this.
TOL_H       = 30.0    # m, height
TOL_XY      = 60.0    # m, ground plan — of the order of a terrain grid cell,
                      # so a reconstructed track stays on the right side of a
                      # ridge instead of cutting through it
COORD_DECIMALS  = 5       # decimals of the position (5 = about 1 m, as in kk7)
# wind measurement (taken over from thermal_strategy.py)
W_MIN_TURNS = 1.5
W_MIN_DUR  = 30     # s
W_TRIM       = 4      # Fixes an jedem Ende abschneiden
K_TURNS_PIECE = 4   # full turns per circling piece: a long climb is cut into
                    # pieces of ~4 turns, each a wind sample of its own
                    # height — the wind profile through the climb instead of
                    # one mean. Measured on 1628 flights: the drift changes
                    # by 2.2 km/h per 100 m of climb in the median (noise
                    # reference, odd vs even turns: 0.3), and the kink of the
                    # circle-centre path follows that change (r = 0.87) but
                    # is the size of the thermal's own wobble (~45 m), so
                    # kinks cannot be located reliably in a single climb;
                    # fixed pieces of 4 turns keep the drift noise of one
                    # piece (~1.4 km/h) well below the signal, and the
                    # profile emerges from many pieces. Runs under 8 turns
                    # stay whole.
W_MAX_DRIFT  = 40.0   # km/h
W_SMOOTH      = 600    # s, smoothing of the wind series between measurements

KINDS = "KG"


# ---------------------------------------------------------------- read IGC
def read_igc(path):
    """-> (day, glider, t0, DataFrame). The glider type comes from the HFGTY
    record of the IGC header; it names no person and is what a sink polar has
    to be fitted per — the per-flight scatter of that quantity is 0.35 m/s,
    far too noisy, while per glider type it settles to about 0.03 m/s."""
    t, la, lo, ab, ag = [], [], [], [], []
    day = 0; glider = ""
    with open(path, "r", encoding="latin-1", errors="ignore") as f:
        for z in f:
            if z.startswith("HFGTY") and not glider:
                glider = re.sub(r"[;|]", " ", z.split(":", 1)[-1]).strip()[:40]
                continue
            if z.startswith("HFDTE") and not day:
                d = re.sub(r"[^0-9]", "", z[5:])
                if len(d) >= 6:
                    j = int(d[4:6]); j += 2000 if j < 80 else 1900
                    day = j * 10000 + int(d[2:4]) * 100 + int(d[0:2])
                continue
            if not z.startswith("B") or len(z) < 35:
                continue
            try:
                t.append(int(z[1:3]) * 3600 + int(z[3:5]) * 60 + int(z[5:7]))
                la.append((int(z[7:9]) + int(z[9:14]) / 60000.0) * (1 if z[14] == "N" else -1))
                lo.append((int(z[15:18]) + int(z[18:23]) / 60000.0) * (1 if z[23] == "E" else -1))
                ab.append(int(z[25:30])); ag.append(int(z[30:35]))
            except ValueError:
                continue
    if len(t) < 60:
        raise ValueError("zu wenige gueltige B-Saetze")
    if not day:
        m = re.search(r"(20\d\d)[-_.]?(\d\d)[-_.]?(\d\d)", Path(path).name)
        day = int(m.group(1)) * 10000 + int(m.group(2)) * 100 + int(m.group(3)) if m else 0
    t = np.array(t, float)
    jump = np.concatenate(([0], np.cumsum(np.diff(t) < 0) * 86400.0))  # midnight rollover
    t += jump
    ab = np.array(ab, float); ag = np.array(ag, float)
    alt = ab if ab.std() > 1 else ag
    return day, glider, float(t[0]), pd.DataFrame(
        dict(t=t - t[0], lat=la, lon=lo, alt=alt))


MAX_GAP_S = 60        # recorder gaps longer than this are not interpolated;
                      # the longest contiguous stretch of the flight is kept


def resample_1hz(df):
    """Resample to 1 Hz, smooth the height, derive vario and ground speed.

    A recorder gap longer than MAX_GAP_S cannot be bridged honestly — an
    interpolated line across it is a teleport that shows up as impossible
    airspeed and absurd air movement. The flight is cut at such gaps and the
    longest contiguous stretch is kept."""
    df = df.sort_values("t").drop_duplicates("t")
    # altitude spikes and GPS warm-up (fixes at 0 m before lock): drop every
    # fix more than 500 m away from the local median of its 31 neighbours
    alt = df["alt"].to_numpy()
    med = pd.Series(alt).rolling(31, center=True, min_periods=1).median().to_numpy()
    df = df[np.abs(alt - med) < 500]
    t = df["t"].to_numpy()
    gap = np.flatnonzero(np.diff(t) > MAX_GAP_S)
    if len(gap):
        gr = np.concatenate(([0], gap + 1, [len(t)]))
        i = int(np.argmax(np.diff(gr)))
        df = df.iloc[gr[i]:gr[i + 1]]
    # t keeps its offset from the first RAW fix, so UTC times stay right
    idx = np.arange(int(df["t"].iloc[0]), int(df["t"].iloc[-1]) + 1)
    out = pd.DataFrame(dict(t=idx.astype(float)))
    for s in ("lat", "lon", "alt"):
        out[s] = np.interp(idx, df["t"].to_numpy(), df[s].to_numpy())
    n = len(out)
    fen = min(15, n if n % 2 else n - 1)
    out["alt_g"] = savgol_filter(out["alt"], max(fen, 5), 2) if n > 5 else out["alt"]
    out["vario"] = pd.Series(np.gradient(out["alt_g"].to_numpy())).rolling(
        15, center=True, min_periods=1).mean().to_numpy()
    la, lo = out["lat"].to_numpy(), out["lon"].to_numpy()
    out["x"] = np.radians(lo - lo[0]) * 6371000.0 * np.cos(np.radians(la.mean()))
    out["y"] = np.radians(la - la[0]) * 6371000.0
    d = np.zeros(n)
    d[1:] = np.hypot(np.diff(out["x"]), np.diff(out["y"]))
    out["gs"] = d * 3.6
    return out


def _abstand(x, y):
    """Largest distance of the points from the straight line end to end."""
    dx, dy = x[-1] - x[0], y[-1] - y[0]
    L = np.hypot(dx, dy)
    if L < 1e-9:
        return float(np.hypot(x - x[0], y - y[0]).max())
    return float((np.abs(dy * (x - x[0]) - dx * (y - y[0])) / L).max())


def dp_corners2(t, hl, x, y):
    """Joint Douglas-Peucker over BOTH measures at once: a span is split at
    the point whose deviation, measured relative to its tolerance (air height
    against TOL_H, ground plan against TOL_XY), is worst — until every piece
    keeps both. One pass, no union of corner sets, no merge step: the
    boundaries are a pure function of the geometry."""
    n = len(t)
    if n < 3:
        return list(range(n))
    keep = np.zeros(n, bool); keep[0] = keep[-1] = True

    def dev(u, v, a, b):
        du, dv = u[b] - u[a], v[b] - v[a]
        L = np.hypot(du, dv)
        if L < 1e-9:
            return np.hypot(u[a + 1:b] - u[a], v[a + 1:b] - v[a])
        return np.abs(dv * (u[a + 1:b] - u[a]) - du * (v[a + 1:b] - v[a])) / L

    stack = [(0, n - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        r1 = dev(t, hl, a, b) / TOL_H
        r2 = dev(x, y, a, b) / TOL_XY
        r = np.maximum(r1, r2)
        i = int(np.argmax(r))
        if r[i] > 1.0:
            k = a + 1 + i
            keep[k] = True
            stack.append((a, k)); stack.append((k, b))
    return list(np.flatnonzero(keep))


def dp_corners(x, y, tol):
    """Douglas-Peucker on a curve; returns the indices of the corners.
    A point survives only if a line without it would be off by more than tol."""
    n = len(x)
    if n < 3:
        return list(range(n))
    keep = np.zeros(n, bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        dx, dy = x[b] - x[a], y[b] - y[a]
        length = np.hypot(dx, dy)
        if length < 1e-9:
            dist = np.hypot(x[a + 1:b] - x[a], y[a + 1:b] - y[a])
        else:
            dist = np.abs(dy * (x[a + 1:b] - x[a])
                          - dx * (y[a + 1:b] - y[a])) / length
        i = int(np.argmax(dist))
        if dist[i] > tol:
            k = a + 1 + i
            keep[k] = True
            stack.append((a, k)); stack.append((k, b))
    return list(np.flatnonzero(keep))


def runs_of(maske):
    """[(a,b)] of the contiguous True stretches."""
    m = np.concatenate(([False], maske.astype(bool), [False]))
    k = np.flatnonzero(np.diff(m.astype(np.int8)))
    return list(zip(k[0::2], k[1::2] - 1))


# ---------------------------------------------------------------- circling
def turn_signal(df):
    la, lo = df["lat"].to_numpy(), df["lon"].to_numpy()
    y = np.sin(np.radians(lo[1:] - lo[:-1])) * np.cos(np.radians(la[1:]))
    x = (np.cos(np.radians(la[:-1])) * np.sin(np.radians(la[1:]))
         - np.sin(np.radians(la[:-1])) * np.cos(np.radians(la[1:]))
         * np.cos(np.radians(lo[1:] - lo[:-1])))
    course = (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0
    dh = np.zeros(len(df))
    dh[2:] = (np.diff(course) + 180.0) % 360.0 - 180.0      # degrees per second
    net = pd.Series(dh).rolling(TURN_WINDOW, center=True,
                                  min_periods=5).sum().to_numpy()
    return dh, np.nan_to_num(net)


def drift_sample(df, a, b, dh):
    """Wind sample of one circling piece, after thermal_strategy: (turns, vx, vy)
    or (turns, nan, nan) if the conditions are not met."""
    turns = abs(float(dh[a:b + 1].sum())) / 360.0
    if b - a < W_MIN_DUR or turns < W_MIN_TURNS:
        return turns, np.nan, np.nan
    a2, b2 = a + W_TRIM, b - W_TRIM
    dt = float(b2 - a2)
    if dt < 20:
        return turns, np.nan, np.nan
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    vx = (x[b2] - x[a2]) / dt
    vy = (y[b2] - y[a2]) / dt
    if np.hypot(vx, vy) * 3.6 > W_MAX_DRIFT:
        return turns, np.nan, np.nan
    return turns, vx, vy


def wind_series(df, dh, net):
    """Wind per second from the circling pieces; only for the state rule."""
    n = len(df)
    tp, ex, ey = [], [], []
    for a, b in runs_of(np.abs(net) > TURN_LIMIT):
        _, vx, vy = drift_sample(df, a, b, dh)
        if np.isfinite(vx):
            tp.append(0.5 * (a + b)); ex.append(vx); ey.append(vy)
    i = np.arange(n, dtype=float)
    if len(tp) >= 2:
        wx, wy = np.interp(i, tp, ex), np.interp(i, tp, ey)
    elif tp:
        wx, wy = np.full(n, ex[0]), np.full(n, ey[0])
    else:
        wx, wy = np.zeros(n), np.zeros(n)
    wx = pd.Series(wx).rolling(W_SMOOTH, center=True, min_periods=1).mean().to_numpy()
    wy = pd.Series(wy).rolling(W_SMOOTH, center=True, min_periods=1).mean().to_numpy()
    return wx, wy, np.hypot(wx, wy) * 3.6


# ---------------------------------------------------------------- states
def circling_per_second(df, net, wind_kmh):
    """First axis: is he circling or not."""
    return (np.abs(net) > TURN_LIMIT) | (df["gs"].to_numpy()
                                           <= WIND_FACTOR * wind_kmh)


def _merge(z, min_s):
    """Absorb runs that are too short into the longer neighbour."""
    z = np.asarray(z).copy()
    while True:
        gr = np.concatenate(([0], np.flatnonzero(np.diff(z)) + 1, [len(z)]))
        dur = np.diff(gr)
        if len(dur) <= 1 or dur.min() >= min_s:
            break
        i = int(np.argmin(dur))
        left = dur[i - 1] if i > 0 else -1
        right = dur[i + 1] if i + 1 < len(dur) else -1
        z[gr[i]:gr[i + 1]] = z[gr[i] - 1] if left >= right else z[gr[i + 1]]
    gr = np.concatenate(([0], np.flatnonzero(np.diff(z)) + 1, [len(z)]))
    return [(int(gr[i]), int(gr[i + 1] - 1), z[gr[i]])
            for i in range(len(gr) - 1)]


ALIAS = {"ozon": "ozone", "adv": "advance", "nivuk": "niviuk"}


def normalise(name):
    """Glider name -> lookup key. Lower-case, punctuation stripped, a few
    manufacturer aliases. Digits are never touched and nothing fuzzy happens:
    "zeolite gt" and "zeolite 2 gt" stay different gliders."""
    t = re.sub(r"[^a-z0-9 ]", " ", name.lower())
    t = re.sub(r"\b(paraglider|glider|wing|pg)\b", " ", t)
    return " ".join(ALIAS.get(w, w) for w in t.split())


def load_polars(path):
    """polars.csv (from polarmaker.py) -> {name: (speeds, sinks)}."""
    tab = {}
    for z in open(path, encoding="utf-8"):
        if z.startswith("#") or z.startswith("glider;") or not z.strip():
            continue
        f = z.rstrip("\n").split(";")
        sp, si = [], []
        for i, cell in enumerate(f[3:]):
            if cell:
                sp.append(25.0 + 5.0 * i); si.append(float(cell))
        if len(sp) >= 3:
            tab[f[0]] = (np.array(sp), np.array(si))
    return tab


def own_sink(df, wx, wy, polars, glider):
    """The glider's own sink per second, from its polar at the momentary
    airspeed. -> (own, source) with source 'g' = the glider's own table row,
    'a' = the _general row, 'n' = no table (constant CALM_SINK)."""
    n = len(df)
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    vx = np.zeros(n); vy = np.zeros(n)
    vx[1:] = np.diff(x); vy[1:] = np.diff(y)
    eig = np.hypot(vx - wx, vy - wy) * 3.6
    eig = np.convolve(eig, np.ones(9) / 9.0, "same")
    if not polars:
        return np.full(n, CALM_SINK), "n", eig
    key = normalise(glider)
    if key in polars:
        sp, si = polars[key]; src = "g"
    else:
        sp, si = polars["_general"]; src = "a"
    return np.interp(np.clip(eig, sp[0], sp[-1]), sp, si), src, eig


def turn_pieces(a, b, dh):
    """Cut a circling run [a, b] into pieces of about K_TURNS_PIECE full
    turns each, boundaries on full turns; fewer than 2*K_TURNS_PIECE turns
    -> one piece."""
    kum = np.abs(dh[a:b + 1]).cumsum()
    turns = float(kum[-1]) / 360.0
    n = int(turns // K_TURNS_PIECE)
    if n < 2:
        return [(a, b)]
    schnitt = [a]
    for k in range(1, n):
        j = a + int(np.searchsorted(kum, turns * k / n * 360.0))
        schnitt.append(min(max(j, schnitt[-1] + 1), b))
    schnitt.append(b + 1)
    return [(schnitt[i], schnitt[i + 1] - 1) for i in range(n)]


def to_runs(df, circ, hl=None, min_s=MIN_RUN_S, dh=None):
    """Runs of equal state.

    The track is first split into circling and straight — a thermal does
    not fall apart at every dip in the climb. A long climb is then cut into
    pieces of K_TURNS_PIECE full turns (turn_pieces), each carrying its own
    drift: the wind profile through the climb. Without dh, circling stays
    one piece.

    A straight run is then divided WHERE THE AIR CHANGES: Douglas-Peucker on
    the air-corrected height hl (height minus the glider's own accumulated
    sink) with TOL_H, and on the ground plan with TOL_XY. Cutting on the raw
    height would put cuts where the pilot changed speed and miss changes of
    the air that his speed change masked.

    Both tolerances are enforced JOINTLY: a span is split at the point whose
    deviation, relative to its tolerance, is worst, until every piece keeps
    both. The boundaries are a pure function of the geometry.
    """
    h = df["alt_g"].to_numpy()
    if hl is None:
        hl = h - CALM_SINK * np.arange(len(h))    # constant own sink
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    t = np.arange(len(h), dtype=float)
    out = []
    for a, b, k in _merge(circ.astype(np.int8), min_s):
        if k:                          # circling: pieces of ~4 turns, K
            for a2, b2 in (turn_pieces(a, b, dh) if dh is not None else [(a, b)]):
                out.append((a2, b2, 0))
            continue
        if b - a < 2:
            out.append((a, b, 1))
            continue
        kn = dp_corners2(t[a:b + 1], hl[a:b + 1], x[a:b + 1], y[a:b + 1])
        for i in range(len(kn) - 1):
            p0, q0 = a + kn[i], a + kn[i + 1]
            e = q0 - 1 if i < len(kn) - 2 else b
            out.append((p0, e, 1))
    return out


def segments(path, polars=None):
    """-> (day, glider, polar_source, segments, end, 1 Hz frame)."""
    day, glider, t0, raw = read_igc(path)
    df = resample_1hz(raw)
    t0 = t0 + float(df["t"].iloc[0])     # the kept block may start later
    dh, net = turn_signal(df)
    wx, wy, wkmh = wind_series(df, dh, net)
    circ = circling_per_second(df, net, wkmh)
    own, src, eig = own_sink(df, wx, wy, polars, glider)
    h = df["alt_g"].to_numpy()
    hl = h - np.cumsum(own)
    la, lo = df["lat"].to_numpy(), df["lon"].to_numpy()
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    schritt = np.zeros(len(df))
    schritt[1:] = np.hypot(np.diff(x), np.diff(y))
    out = []
    for a, b, k in to_runs(df, circ, hl, dh=dh):
        kind = KINDS[k]
        s = dict(kind=kind, t=int(round(t0 + a)) % 86400, h=int(round(h[a])),
                 la=float(la[a]), lo=float(lo[a]), a=a, b=b)
        if kind in "Kk":
            kr, vx, vy = drift_sample(df, a, b, dh)
            s["turns"] = round(kr, 1)
            if np.isfinite(vx):
                s["drift"] = round(float(np.hypot(vx, vy)) * 3.6, 1)
                # meteorological: the direction the air comes FROM
                s["dir"] = int(round((np.degrees(np.arctan2(-vx, -vy)) + 360) % 360))
            else:
                s["drift"] = None; s["dir"] = None
        else:
            # the movement of the air over this piece, m/s — own sink at the
            # flown airspeed already removed; ready to colour a map with
            w = round((hl[b] - hl[a]) / max(b - a, 1), 1)
            s["w"] = 0.0 if w == 0 else w      # never "-0.0" in the line
            # true airspeed (mean of the seconds — NOT the chord over the
            # duration, which a curved path underestimates) ...
            s["v"] = int(round(eig[a + 1:b + 1].mean())) if b > a else 0
            # ... and how much longer the flown path was than the chord, in %
            sehne = float(np.hypot(x[b] - x[a], y[b] - y[a]))
            weg = float(schritt[a + 1:b + 1].sum())
            s["z"] = int(min(round(100 * (weg / sehne - 1)), 999)) \
                if sehne > 1 else 0
            s["z"] = max(s["z"], 0)
        out.append(s)
    end = (int(round(t0 + len(df) - 1)) % 86400, int(round(h[-1])),
            float(la[-1]), float(lo[-1]))
    return day, glider, src, out, end, df


def line(kennung, glider, pol, day, segs, end):
    st = []
    for s in segs:
        Q = 10**COORD_DECIMALS      # print EXACTLY what line_delta rounds to
        b = (f"{s['kind']},{s['t']},{s['h']},"
             f"{round(s['la']*Q)/Q:.{COORD_DECIMALS}f},"
             f"{round(s['lo']*Q)/Q:.{COORD_DECIMALS}f}")
        if s["kind"] in "Kk":
            d = "" if s["drift"] is None else f"{s['drift']:.1f}"
            r = "" if s["dir"] is None else str(s["dir"])
            b += f",{s['turns']:.1f},{d},{r}"
        else:
            b += f",{s['w']:.1f},{s['v']},{s['z']}"
        st.append(b)
    return (f"{kennung};{glider};{pol};{day};" + "|".join(st)
            + f";{end[0]},{end[1]},"
              f"{round(end[2]*Q)/Q:.{COORD_DECIMALS}f},"
              f"{round(end[3]*Q)/Q:.{COORD_DECIMALS}f}")


def line_delta(kennung, glider, pol, day, segs, end):
    """The same line, but every number as the difference to the previous one.
    Nothing is lost — read_line_delta() gives back exactly what read_line()
    gives — the numbers merely get short."""
    st, v = [], None
    for s_ in segs:
        la, lo = round(s_["la"] * 10**COORD_DECIMALS), round(s_["lo"] * 10**COORD_DECIMALS)
        if v is None:
            b = f"{s_['kind']}{s_['t']},{s_['h']},{la},{lo}"
        else:
            b = (f"{s_['kind']}{s_['t']-v[0]},{s_['h']-v[1]},"
                 f"{la-v[2]},{lo-v[3]}")
        if s_["kind"] in "Kk":
            d = "" if s_["drift"] is None else f"{s_['drift']:.1f}"
            r = "" if s_["dir"] is None else str(s_["dir"])
            b += f",{s_['turns']:.1f},{d},{r}"
        else:
            b += f",{s_['w']:.1f},{s_['v']},{s_['z']}"
        st.append(b); v = (s_["t"], s_["h"], la, lo)
    e = (end[0] - v[0], end[1] - v[1],
         round(end[2] * 10**COORD_DECIMALS) - v[2], round(end[3] * 10**COORD_DECIMALS) - v[3])
    return f"{kennung};{glider};{pol};{day};" + "|".join(st) + ";" + ",".join(str(z) for z in e)


def read_line_delta(zl):
    """Counterpart of line_delta()."""
    kennung, glider, pol, day, sr, er = zl.rstrip("\n").split(";")
    segs, v = [], None
    for stk in sr.split("|"):
        f = stk.split(",")
        t, h, la, lo = int(f[0][1:]), int(f[1]), int(f[2]), int(f[3])
        if v is not None:
            t += v[0]; h += v[1]; la += v[2]; lo += v[3]
        v = (t, h, la, lo)
        s_ = dict(kind=f[0][0], t=t, h=h, la=la / 10**COORD_DECIMALS, lo=lo / 10**COORD_DECIMALS)
        if s_["kind"] in "Kk":
            s_["turns"] = float(f[4])
            s_["drift"] = float(f[5]) if f[5] else None
            s_["dir"] = int(f[6]) if f[6] else None
        else:
            s_["w"] = float(f[4]); s_["v"] = int(f[5]); s_["z"] = int(f[6])
        segs.append(s_)
    f = [int(z) for z in er.split(",")]
    end = (v[0] + f[0], v[1] + f[1], (v[2] + f[2]) / 10**COORD_DECIMALS, (v[3] + f[3]) / 10**COORD_DECIMALS)
    return kennung, glider, pol, int(day), segs, end


def read_line(zl):
    """Counterpart of line(): line -> (id, glider, pol, day, segments, end)."""
    kennung, glider, pol, day, sr, er = zl.rstrip("\n").split(";")
    segs = []
    for st in sr.split("|"):
        f = st.split(",")
        s = dict(kind=f[0], t=int(f[1]), h=int(f[2]), la=float(f[3]), lo=float(f[4]))
        if s["kind"] in "Kk":
            s["turns"] = float(f[5])
            s["drift"] = float(f[6]) if f[6] else None
            s["dir"] = int(f[7]) if f[7] else None
        else:
            s["w"] = float(f[5]); s["v"] = int(f[6]); s["z"] = int(f[7])
        segs.append(s)
    f = er.split(",")
    return kennung, glider, pol, int(day), segs, (int(f[0]), int(f[1]), float(f[2]), float(f[3]))


if __name__ == "__main__":
    import hashlib
    argv = sys.argv[1:]
    short = "--delta" in argv
    if short:
        argv.remove("--delta")
    part = [a for a in argv if a.startswith("--part")]
    k, n = 1, 1
    if part:
        argv.remove(part[0])
        k, n = (int(x) for x in part[0].split("=")[1].split("/"))
    kennung = "P"
    if "--id" in argv:
        i = argv.index("--id"); kennung = argv[i + 1]; del argv[i:i + 2]
    poldat = None
    if "--polars" in argv:
        i = argv.index("--polars"); poldat = argv[i + 1]; del argv[i:i + 2]
    elif (Path(__file__).resolve().parent / "polars.csv").exists():
        poldat = str(Path(__file__).resolve().parent / "polars.csv")
    polars = load_polars(poldat) if poldat else None
    if poldat:
        pruef = hashlib.sha1(open(poldat, "rb").read()).hexdigest()[:12]
        print(f"# flightstates 2.0 polars={Path(poldat).name} sha1={pruef}")
    else:
        print("# flightstates 2.0 polars=none (constant own sink "
              f"{CALM_SINK} m/s)")
    dateien = []
    for a in argv:
        trf = _glob.glob(a, recursive=True)
        dateien += trf if trf else [a]
    einmal = {}
    for d in sorted(dateien):
        einmal.setdefault(Path(d).name, d)
    dateien = sorted(einmal.values())[k - 1::n]
    for p in dateien:
        try:
            day, glider, src, segs, end, _ = segments(p, polars)
            k_ = Path(p).stem if kennung == "name" else kennung
            f = line_delta if short else line
            print(f(k_, glider, src, day, segs, end))
        except Exception as e:
            print(f"# {Path(p).name}: {type(e).__name__}: {e}", file=sys.stderr)
