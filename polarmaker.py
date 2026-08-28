#!/usr/bin/env python3
"""
polarmaker.py — run 1 of 2: measure the effective sink polar per glider type.

Reads IGC files and writes ONE small table (polars.csv) and nothing else:
per glider type, the glider's sink at airspeeds 25-65 km/h, measured from the
glide legs of all its flights. No positions, no times, no pilots, no flight
data of any kind leave this script — the table is a few kilobytes and can be
inspected line by line before anything else is run.

Why this table exists: run 2 (flightstates.py) cuts a flight into segments
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

    python3 polarmaker.py "flights/**/*.IGC" --out polars.csv
    python3 polarmaker.py "flights/**/*.IGC" --part 2/4 --out part2.csv
    python3 polarmaker.py --join "part*.csv" --out polars.csv

Dependencies: numpy, pandas, scipy, and flightstates.py beside this file.
"""
import argparse, glob, re, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import flightstates as fs                                    # noqa: E402

SPEEDS = list(range(25, 66, 5))          # band centres, km/h (width 5)
MIN_SEC_CELL = 3000                      # seconds needed before a cell is written

ALIAS = {"ozon": "ozone", "adv": "advance", "nivuk": "niviuk"}


def normalise(name):
    """Lower-case, strip punctuation, alias manufacturers. Digits untouched."""
    s = re.sub(r"[^a-z0-9 ]", " ", name.lower())
    s = re.sub(r"\b(paraglider|glider|wing|pg)\b", " ", s)
    return " ".join(ALIAS.get(w, w) for w in s.split())


def seconds_of(path):
    """(glider, airspeeds, varios) — one pair per non-circling second."""
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
    m = (~circ) & (eig > 18) & (eig < 75) & (v > -6) & (v < 4)
    return normalise(glider), eig[m].astype(np.float32), v[m].astype(np.float32)


def mode_of(v, width=0.05):
    """The most common value, parabola-refined over a smoothed histogram."""
    k = np.arange(-6, 4.001, width)
    h, _ = np.histogram(v, bins=k)
    h = np.convolve(h, np.ones(5) / 5.0, "same")
    i = int(np.argmax(h))
    if 0 < i < len(h) - 1:
        y0, y1, y2 = h[i - 1], h[i], h[i + 1]
        d = y0 - 2 * y1 + y2
        c = 0.5 * (y0 - y2) / d if abs(d) > 1e-9 else 0.0
    else:
        c = 0.0
    return k[i] + width * (0.5 + c)


def write_table(dest, sammel, min_flights):
    rows = []
    for name, (fl, ee, vv) in sammel.items():
        if not name or len(fl) < min_flights or not len(ee):
            continue
        e = np.concatenate(ee); v = np.concatenate(vv)
        werte = []
        for c in SPEEDS:
            m = (e >= c - 2.5) & (e < c + 2.5)
            werte.append(mode_of(v[m]) if m.sum() >= MIN_SEC_CELL else None)
        # Below 35 km/h the most common air is ridge lift, not still air, and
        # the mode can land on the lift instead of the glider. A slow cell is
        # kept only if it agrees with the 35 band to 0.25 m/s.
        anker = werte[SPEEDS.index(35)]
        for i, c in enumerate(SPEEDS):
            if c < 35 and werte[i] is not None:
                if anker is None or werte[i] > anker + 0.25:
                    werte[i] = None
        cells = ["" if w is None else f"{w:.2f}" for w in werte]
        if sum(1 for c in cells if c) >= 3:
            rows.append((name, len(fl), int(len(e)), cells))
    rows.sort(key=lambda r: -r[1])
    with open(dest, "w", encoding="utf-8") as f:
        f.write("# flightstates polar table 1.1\n")
        f.write("# own sink [m/s] of the glider in near-still air (mode per "
                "airspeed band [km/h])\n")
        f.write("glider;flights;seconds;" + ";".join(str(c) for c in SPEEDS) + "\n")
        for name, nf, nl, cells in rows:
            f.write(f"{name};{nf};{nl};" + ";".join(cells) + "\n")
    print(f"{dest}: {len(rows)} rows, {Path(dest).stat().st_size} bytes")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("globs", nargs="*")
    ap.add_argument("--out", default="polars.csv")
    ap.add_argument("--part", default="1/1", help="k/n for parallel runs")
    ap.add_argument("--min", type=int, default=20,
                    help="flights needed before a glider gets a row (use 50+ "
                         "on a large archive)")
    ap.add_argument("--join", default=None,
                    help="glob of part tables to merge instead of reading IGCs")
    a = ap.parse_args()

    if a.join:
        # merge part tables: per cell the second-weighted mean of the modes
        agg = {}
        for part in sorted(glob.glob(a.join)):
            for z in open(part, encoding="utf-8"):
                if z.startswith("#") or z.startswith("glider;"):
                    continue
                f = z.rstrip("\n").split(";")
                name, nf, ns = f[0], int(f[1]), int(f[2])
                e = agg.setdefault(name, [0, 0, [[0.0, 0] for _ in SPEEDS]])
                e[0] += nf; e[1] += ns
                for i, cell in enumerate(f[3:]):
                    if cell:
                        e[2][i][0] += float(cell) * ns; e[2][i][1] += ns
        with open(a.out, "w", encoding="utf-8") as f:
            f.write("# flightstates polar table 1.1 (joined)\n")
            f.write("glider;flights;seconds;" + ";".join(str(c) for c in SPEEDS) + "\n")
            for name, (nf, ns, cells) in sorted(agg.items(), key=lambda z: -z[1][0]):
                if nf < a.min:
                    continue
                f.write(f"{name};{nf};{ns};" + ";".join(
                    f"{c[0]/c[1]:.2f}" if c[1] else "" for c in cells) + "\n")
        print(f"{a.out} written")
        return

    files = []
    for g in a.globs:
        files += glob.glob(g, recursive=True)
    einmal = {}
    for d in sorted(files):
        einmal.setdefault(Path(d).name, d)
    files = sorted(einmal.values())
    k, n = (int(x) for x in a.part.split("/"))
    files = files[k - 1::n]

    sammel = {}                # name -> (set of flight idx, [eig], [vario])
    ae, av = [], []
    for i, p in enumerate(files, 1):
        try:
            name, e, v = seconds_of(p)
        except Exception as err:
            print(f"# {Path(p).name}: {type(err).__name__}: {err}", file=sys.stderr)
            continue
        fl, ee, vv = sammel.setdefault(name, (set(), [], []))
        fl.add(i); ee.append(e); vv.append(v)
        ae.append(e); av.append(v)
        if i % 100 == 0:
            print(f"  {i}/{len(files)}", file=sys.stderr)
    sammel["_general"] = (set(range(len(files))), ae, av)
    write_table(a.out, sammel, a.min)


if __name__ == "__main__":
    main()
