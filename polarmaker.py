#!/usr/bin/env python3
"""
polarmaker.py — measure the effective sink polar per glider type.

Reads IGC files and writes three files: polars.csv (the table — per glider
type, the glider's sink at airspeeds 25-65 km/h, measured from the glide
legs of all its flights), polars_stat.csv (the statistics behind every
value) and polars.npz (the raw per-flight histograms). No positions, no
times, no pilots, no flight data of any kind leave this script — the text
files can be inspected line by line, the .npz holds nothing but count
tables per glider type.

Why the table exists: the segmenter (flightstates.py) cuts a flight into segments
where the AIR changes, not where the vario changes. For that it must subtract
the glider's own sink at the momentary airspeed — and that curve differs from
glider to glider. Estimated from a single flight it scatters by 0.21 m/s;
averaged over 50+ flights of the same type it settles to +-0.03.

Method: circling episodes give the wind (drift over full turns). Every
non-circling SECOND contributes one pair (airspeed, vario); airspeed is the
ground velocity minus the wind vector. Pairs are pooled per NORMALISED glider
name, and per 5 km/h airspeed band the MODE of the vario is taken — the most
common value, not the median. The mode matters: the most common air on a
glide is near-still air, so the mode strips the lift and sink out and lands
on the glider's own sink, while a median would keep the lift in and make
every polar look a shade too good. Checked against two independent anchors
on 774 flights, the mode reproduces them to a few hundredths.

Names: lower-cased, punctuation stripped, a few manufacturer aliases. Digits
are never touched and no fuzzy matching is done — "zeolite gt" and
"zeolite 2 gt" are different gliders. Typos stay separate rows and fall
below the --min threshold on their own.

Discarding spoiled cells: two rules keep venue-biased corpora (e.g. a
competition wing flown only at one ridge site, in lift, on a task clock)
from writing air into the polar. (1) Physics floor SINK_MIN: no paraglider
sinks less than 0.85 m/s in straight flight — a cell above that measured
lift, not the glider. (2) Bimodality gates: where the bootstrap whisker is
wider than BREITE_MAX, or more than BODEN_MAX of the bootstrap draws hits
the physics floor, the mode jumps between a still-air and a lift peak and
the cell is discarded as unstable. Whole-flight filtering was measured and rejected: the share of
rising slow seconds overlaps too much between clean and spoiled corpora
(0.31/0.38/0.44 for zeolite2/zeno2/enzo3), and it would discard the
flights' perfectly good fast-band seconds with them.

Whiskers: polars_stat.csv records EVERY band of EVERY glider — kept or
discarded, with the reason — and the 90 % interval of a flight-cluster
bootstrap (see bootstrap_stats). polars.csv itself stays byte-identical in
format; nothing downstream changes.

Pots in the .npz (since 1.2): the mode assumes the most common air on a
glide is near-still air. That assumption is weakest at the fast end, where
pilots push the speed bar BECAUSE the air sinks. To make it checkable
later, every non-circling second is sorted into one of three pots before it
is counted, and the .npz keeps a histogram per pot and flight: h = all
seconds (what the table is built from — unchanged), q = QUIET (steady
vario: std over QUIET_WIN_S below QUIET_STD_MS — still air does not
fluctuate, thermic air does; the level of the vario is never looked at,
only its steadiness — no circling within QUIET_NOCIRC_S, not right after
a bar change), a = ACCEL (the ACCEL_AFTER_S after the airspeed, 15-s mean,
changed by ACCEL_KMH within ACCEL_WIN_S). Rest = h - q - a. The table and
the whisker file are computed from h exactly as before; q and a are there
to be looked at. Counts only, as before — no positions, no times, no
pilots leave the script.
Parts dump per-flight histograms (.npz); --join merges them and builds the
table with the full statistics and gates — identical machinery to a single
pass. Memory stays small either way: per flight only a 4 KB histogram is
kept, never the raw seconds.

    python3 polarmaker.py "flights/**/*.IGC" --out polars.csv
    python3 polarmaker.py "flights/**/*.IGC" --part 2/4 --out teil2.npz
    python3 polarmaker.py --join "teil*.npz" --out polars.csv

Dependencies: numpy, pandas, scipy, and flightstates.py beside this file.
"""
import argparse, glob, json, re, sys, zlib
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import flightstates as fs                                    # noqa: E402

SPEEDS = list(range(25, 66, 5))          # band centres, km/h (width 5)
MIN_SEC_CELL = 3000                      # seconds needed before a cell is written
SINK_MIN = -0.85                         # physics floor: no paraglider sinks less
                                         # in straight flight — a cell above this
                                         # measured the AIR (ridge lift), not the
                                         # glider, and is discarded
BREITE_MAX = 0.75                        # bimodality gate: a cell whose 90 %
                                         # bootstrap interval is wider than this
                                         # has a mode that jumps between two
                                         # peaks (still air vs lift) and is
                                         # discarded. Measured gap in the
                                         # archive: healthy cells <= 0.68,
                                         # pathological ones >= 0.85
BODEN_MAX = 0.05                         # stability gate: if more than this
                                         # share of bootstrap draws hits the
                                         # physics floor, the still-air peak is
                                         # not dominant enough to trust the
                                         # cell. Measured gap: healthy cells
                                         # 0.00, spoiled ones 0.07-0.33

QUIET_NOCIRC_S = 120                     # quiet pot: no circling within +-N s
QUIET_STD_MS = 0.35                      # ... vario std over QUIET_WIN_S below this
QUIET_WIN_S = 61                         # ... (odd, centred)
ACCEL_KMH = 10.0                         # accel pot: airspeed (15-s mean) changed by >= this
ACCEL_WIN_S = 30                         # ... over this many seconds. GPS airspeed
                                         # jitters 3-4 km/h second to second; anything
                                         # finer than 15-s means flags noise, not the bar
ACCEL_AFTER_S = 30                       # ... marks the N s after the change
POTS = ("all", "quiet", "accel")

ALIAS = {"ozon": "ozone", "adv": "advance", "nivuk": "niviuk"}


def normalise(name):
    """Lower-case, strip punctuation, alias manufacturers. Digits untouched."""
    s = re.sub(r"[^a-z0-9 ]", " ", name.lower())
    s = re.sub(r"\b(paraglider|glider|wing|pg)\b", " ", s)
    return " ".join(ALIAS.get(w, w) for w in s.split())


def pots_of(eig, circ, vario):
    """Pot of every second: 0 rest, 1 quiet, 2 accel (see module docstring).
    A label per second, nothing is stored."""
    import pandas as pd
    n = len(eig)
    e15 = np.convolve(eig, np.ones(15) / 15.0, "same")
    d = np.zeros(n)
    d[ACCEL_WIN_S:] = np.abs(e15[ACCEL_WIN_S:] - e15[:-ACCEL_WIN_S])
    accel = np.zeros(n, bool)
    for i in np.flatnonzero(d >= ACCEL_KMH):
        accel[i:i + ACCEL_AFTER_S + 1] = True
    ci = np.flatnonzero(circ)
    idx = np.arange(n)
    if len(ci):
        pos = np.searchsorted(ci, idx)
        prev = np.where(pos > 0, ci[np.clip(pos - 1, 0, len(ci) - 1)], -10 ** 9)
        nxt = np.where(pos < len(ci), ci[np.clip(pos, 0, len(ci) - 1)], 10 ** 9)
        abstand = np.minimum(idx - prev, nxt - idx)
    else:
        abstand = np.full(n, 10 ** 9)
    std = pd.Series(vario).rolling(QUIET_WIN_S, center=True, min_periods=QUIET_WIN_S // 2).std().to_numpy()
    quiet = (abstand > QUIET_NOCIRC_S) & (std < QUIET_STD_MS) & ~accel
    pot = np.zeros(n, np.int8)
    pot[accel] = 2
    pot[quiet] = 1
    return pot


def seconds_of(path):
    """(glider, airspeeds, varios, pots) — one triple per non-circling second."""
    day, glider, t0, raw = fs.read_igc(path)
    df = fs.resample_1hz(raw)
    dh, net = fs.turn_signal(df)
    wx, wy, wk = fs.wind_series(df, dh, net)
    circ = fs.circling_per_second(df, net, wk)
    n = len(df)
    x = df["x"].to_numpy(); y = df["y"].to_numpy()
    vx = np.zeros(n); vy = np.zeros(n)
    vx[1:] = np.diff(x); vy[1:] = np.diff(y)
    eig = np.hypot(vx - wx, vy - wy) * 3.6
    v = df["vario"].to_numpy()
    pot = pots_of(eig, circ, v)
    m = (~circ) & (eig > 18) & (eig < 75) & (v > -6) & (v < 4)
    return normalise(glider), eig[m].astype(np.float32), v[m].astype(np.float32), pot[m]


KANTEN = np.arange(-6, 4.001, 0.05)      # the one histogram grid everything uses


def mode_of(v, width=0.05):
    """The most common value, parabola-refined over a smoothed histogram."""
    h, _ = np.histogram(v, bins=KANTEN)
    return mode_from_hist(h, KANTEN, width)


def hist_of(e, v, pot=None):
    """One flight -> (pot x band x bin) uint16 histograms: [0] all seconds,
    [1] quiet, [2] accel. All later maths — mode, bootstrap, gates — needs
    only [0]; the raw seconds are never kept, so memory stays a few KB per
    flight even on a 100k-flight archive. uint16 is safe: no flight has
    65535 seconds in one band and bin."""
    H = np.zeros((len(POTS), len(SPEEDS), len(KANTEN) - 1), np.uint16)
    if pot is None:
        pot = np.zeros(len(e), np.int8)
    for i, c in enumerate(SPEEDS):
        m = (e >= c - 2.5) & (e < c + 2.5)
        if m.any():
            H[0, i] = np.histogram(v[m], bins=KANTEN)[0]
            H[1, i] = np.histogram(v[m & (pot == 1)], bins=KANTEN)[0]
            H[2, i] = np.histogram(v[m & (pot == 2)], bins=KANTEN)[0]
    return H


def mode_from_hist(h, k, width=0.05):
    """mode_of, but starting from a ready-made histogram on the grid k."""
    h = np.convolve(h, np.ones(5) / 5.0, "same")
    i = int(np.argmax(h))
    if 0 < i < len(h) - 1:
        y0, y1, y2 = h[i - 1], h[i], h[i + 1]
        d = y0 - 2 * y1 + y2
        c = 0.5 * (y0 - y2) / d if abs(d) > 1e-9 else 0.0
    else:
        c = 0.0
    return k[i] + width * (0.5 + c)


B_BOOT = 200                             # bootstrap draws for the whiskers


def bootstrap_stats(HL, rng):
    """Flight-cluster bootstrap of the per-band mode -> whisker data.

    Seconds within one flight share the same air, trim and speed-bar style
    and are heavily correlated; resampling seconds would give absurdly
    narrow intervals. The resampling unit is therefore the FLIGHT: draw
    flights with replacement, sum their per-band histograms, take the mode
    again — with the same MIN_SEC_CELL and slow-band anchor rule applied
    per draw. Returns per band (lo5, hi95, kept_share) over B_BOOT draws,
    or None where the cell survived too few draws to be quoted.
    """
    width = 0.05
    k = KANTEN
    H = np.stack(HL)                     # (flights x bands x bins), uint16
    F = H.shape[0]
    nb = len(SPEEDS)
    zieh = [[] for _ in SPEEDS]
    boden = [0] * nb                     # draws in which the mode hit the floor
    for _ in range(B_BOOT):
        w = np.bincount(rng.integers(0, F, F), minlength=F).astype(np.int64)
        hs = np.tensordot(w, H, axes=1)  # = H[idx].sum(0), but fast at 100k
        werte = []
        for i in range(nb):
            werte.append(mode_from_hist(hs[i], k, width)
                         if hs[i].sum() >= MIN_SEC_CELL else None)
        for i in range(nb):
            if werte[i] is not None and werte[i] > SINK_MIN:
                werte[i] = None                       # physics floor, per draw
                boden[i] += 1
        anker = werte[SPEEDS.index(35)]
        for i, c in enumerate(SPEEDS):
            if c < 35 and werte[i] is not None:
                if anker is None or werte[i] > anker + 0.25:
                    werte[i] = None
        for i in range(nb):
            if werte[i] is not None:
                zieh[i].append(werte[i])
    aus = []
    for i in range(nb):
        z = sorted(zieh[i])
        if len(z) >= B_BOOT // 5:
            aus.append((z[int(0.05 * (len(z) - 1))],
                        z[int(0.95 * (len(z) - 1))],
                        len(z) / B_BOOT,
                        boden[i] / B_BOOT))
        else:
            aus.append(None)
    return aus


def dump_schreiben(dest, sammel):
    """Per-flight histograms of every glider -> one compressed .npz.
    Counts only: no positions, no times, no pilots."""
    arrs = {"namen": np.array(json.dumps(list(sammel))), "toepfe": np.array(json.dumps(POTS))}
    for i, k2 in enumerate(sammel):
        nf, HL, LE = sammel[k2]
        H = np.stack(HL)                                 # (flights x pots x bands x bins)
        arrs[f"h{i}"] = H[:, 0]                          # all seconds — same as before
        arrs[f"q{i}"] = H[:, 1]                          # quiet pot
        arrs[f"a{i}"] = H[:, 2]                          # accel pot
        arrs[f"l{i}"] = np.array(LE, np.int64)
    np.savez_compressed(dest, **arrs)
    print(f"{dest}: raw per-flight histograms, {Path(dest).stat().st_size} bytes")


def write_table(dest, sammel, min_flights):
    rows = []
    stat = []                        # whisker sidecar, one line per kept cell
    for name, (nf, HL3, LE) in sammel.items():
        if not name or nf < min_flights or not HL3:
            continue
        HL = [np.asarray(h)[0] for h in HL3]            # all seconds — the table as before
        Hs = np.zeros((len(SPEEDS), len(KANTEN) - 1), np.int64)
        for h in HL:
            Hs += h
        werte = []
        for i, c in enumerate(SPEEDS):
            werte.append(mode_from_hist(Hs[i], KANTEN)
                         if Hs[i].sum() >= MIN_SEC_CELL else None)
        roh = list(werte)          # before any gate — the sidecar keeps all
        grund = ["" if w is not None else "duenn" for w in werte]
        # Physics floor: a mode above SINK_MIN measured the air, not the
        # glider (competition flights ridge-soaring at slow speed do this).
        for i in range(len(SPEEDS)):
            if werte[i] is not None and werte[i] > SINK_MIN:
                werte[i] = None
                grund[i] = "boden"
        # Below 35 km/h the most common air is ridge lift, not still air, and
        # the mode can land on the lift instead of the glider. A slow cell is
        # kept only if it agrees with the 35 band to 0.25 m/s.
        anker = werte[SPEEDS.index(35)]
        for i, c in enumerate(SPEEDS):
            if c < 35 and werte[i] is not None:
                if anker is None or werte[i] > anker + 0.25:
                    werte[i] = None
                    grund[i] = "anker"
        # Bimodality gates: where the bootstrap interval is wider than
        # BREITE_MAX, or more than BODEN_MAX of the draws hit the physics
        # floor, the mode jumps between a still-air and a lift peak — the
        # cell is not a property of the glider and is discarded.
        # the same corpus gives the same whiskers, no matter how it was
        # processed: seed from the glider name, and put the flights into a
        # canonical order first (a joined run lists them differently than a
        # single pass)
        ordnung = sorted(range(len(HL)),
                         key=lambda f: (LE[f], HL[f].tobytes()))
        HL = [HL[f] for f in ordnung]
        rng = np.random.default_rng(zlib.crc32(name.encode()))
        st = bootstrap_stats(HL, rng)
        for i in range(len(SPEEDS)):
            if werte[i] is not None and st[i] is not None:
                lo, hi, _teil, geboden = st[i]
                if hi - lo > BREITE_MAX or geboden > BODEN_MAX:
                    werte[i] = None
                    grund[i] = "breite" if hi - lo > BREITE_MAX else "instabil"
        cells = ["" if w is None else f"{w:.2f}" for w in werte]
        if sum(1 for c in cells if c) >= 3:
            rows.append((name, nf, int(sum(LE)), cells))
        # The sidecar keeps EVERYTHING — every band of every glider, kept or
        # discarded, with the reason. The final pick/blacklist of glider
        # types is decided later, from this file (or the .npz dumps), with
        # all numbers on the table; the main table is only today's default.
        for i, c in enumerate(SPEEDS):
            m0 = "" if roh[i] is None else f"{roh[i]:.2f}"
            if st[i] is not None:
                lo, hi, teil, geboden = st[i]
                rest = f"{lo:.2f};{hi:.2f};{teil:.2f};{geboden:.2f}"
            else:
                rest = ";;;"
            stat.append(f"{name};{nf};{c};{m0};{rest};"
                        f"{int(Hs[i].sum())};{grund[i] or 'ok'}")
    rows.sort(key=lambda r: -r[1])
    with open(dest, "w", encoding="utf-8") as f:
        f.write("# flightstates polar table 1.1\n")
        f.write("# own sink [m/s] of the glider in near-still air (mode per "
                "airspeed band [km/h])\n")
        f.write("glider;flights;seconds;" + ";".join(str(c) for c in SPEEDS) + "\n")
        for name, nf, nl, cells in rows:
            f.write(f"{name};{nf};{nl};" + ";".join(cells) + "\n")
    print(f"{dest}: {len(rows)} rows, {Path(dest).stat().st_size} bytes")
    # Sidecar: the COMPLETE statistical record — every band of every glider,
    # kept or discarded (status ok/boden/anker/breite/instabil/duenn), with
    # the 90 % flight-bootstrap interval. A later pick or blacklist of
    # glider types is decided from this file; nothing is thrown away here.
    stat_dest = re.sub(r"\.csv$", "", dest) + "_stat.csv"
    with open(stat_dest, "w", encoding="utf-8") as f:
        f.write("# flightstates polar whiskers 1.0 — flight-cluster bootstrap, "
                f"{B_BOOT} draws, 90 % interval\n")
        f.write("glider;flights;band_kmh;mode;lo5;hi95;kept;floored;"
                "seconds_band;status\n")
        for z in stat:
            f.write(z + "\n")
    print(f"{stat_dest}: {len(stat)} cells, {Path(stat_dest).stat().st_size} bytes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("globs", nargs="*")
    ap.add_argument("--out", default="polars.csv")
    ap.add_argument("--part", default="1/1", help="k/n for parallel runs")
    ap.add_argument("--min", type=int, default=20,
                    help="pre-filter: flights needed before a glider gets a "
                         "row. Leave it — the statistics files always cover "
                         "every glider, whatever this is set to")
    ap.add_argument("--join", default=None,
                    help="glob of part tables to merge instead of reading IGCs")
    a = ap.parse_args()

    if a.join:
        # merge part dumps (.npz with per-flight histograms) and build the
        # table with the FULL machinery — statistics and gates identical to
        # a single pass. (Old CSV part tables cannot be merged this way.)
        sammel = {}
        for part in sorted(glob.glob(a.join)):
            z = np.load(part)
            keys = json.loads(str(z["namen"]))
            for i, k in enumerate(keys):
                H = z[f"h{i}"]; L = z[f"l{i}"]
                Q = z[f"q{i}"] if f"q{i}" in z else np.zeros_like(H)   # dumps of 1.1: no pots
                A = z[f"a{i}"] if f"a{i}" in z else np.zeros_like(H)
                g = sammel.setdefault(k, [0, [], []])
                g[0] += H.shape[0]
                g[1].extend(np.stack([H[j], Q[j], A[j]]) for j in range(H.shape[0]))
                g[2].extend(int(x) for x in L)
        write_table(a.out, sammel, a.min)
        # the merged raw record, same as a single pass would ship — but
        # never overwrite one of the input dumps (e.g. --join probe.npz)
        dnpz = re.sub(r"\.csv$", "", a.out) + ".npz"
        import os
        if os.path.abspath(dnpz) not in {os.path.abspath(x)
                                         for x in glob.glob(a.join)}:
            dump_schreiben(dnpz, sammel)
        return

    files = []
    for g in a.globs:
        files += glob.glob(g, recursive=True)
    # duplicates: same NAME AND SIZE is the same file twice. (Name alone
    # would silently drop different flights that happen to share a file
    # name — real in large foreign archives.)
    einmal = {}
    for d in sorted(files):
        try:
            gr = Path(d).stat().st_size
        except OSError:
            continue
        einmal.setdefault((Path(d).name, gr), d)
    files = sorted(einmal.values())
    k, n = (int(x) for x in a.part.split("/"))
    files = files[k - 1::n]

    sammel = {}                # name -> [flights, [flight histograms], [seconds]]
    for i, p in enumerate(files, 1):
        try:
            name, e, v, pot = seconds_of(p)
        except Exception as err:
            print(f"# {Path(p).name}: {type(err).__name__}: {err}", file=sys.stderr)
            continue
        H = hist_of(e, v, pot)
        for nm in (name, "_general"):
            g = sammel.setdefault(nm, [0, [], []])
            g[0] += 1
            g[1].append(H)             # same object twice — no copy
            g[2].append(int(len(e)))
        if i % 100 == 0:
            print(f"  {i}/{len(files)}", file=sys.stderr)
    if n > 1:
        # a part writes no table — it dumps its per-flight histograms, and
        # --join builds the table with full statistics and gates.
        dest = a.out if a.out.endswith(".npz") else re.sub(r"\.csv$", "", a.out) + ".npz"
        dump_schreiben(dest, sammel)
        print(f"{dest}: part {k}/{n}, {len(files)} files — merge with --join")
        return
    write_table(a.out, sammel, a.min)
    # The raw statistical record always ships with the table: per-flight
    # histograms of every glider, below --min included. From this file every
    # later selection or re-check can be computed without the IGCs.
    dump_schreiben(re.sub(r"\.csv$", "", a.out) + ".npz", sammel)


if __name__ == "__main__":
    main()
