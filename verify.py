"""verify.py — the same picture twice: on top from the real IGC file, below
reconstructed from the text line alone.

  python3 verify.py flight.IGC --out check.png

Both blocks show the same five rows with the same colours and the same
definitions. The upper block computes from all 1 Hz points of the track, the
lower one knows only the text line — per segment the kind, the second, the
height, the position, and for circling segments the turn count and the drift.
That is why the lower block is angular: between two segment starts it joins
the points by a straight line.

It also prints the numeric comparison: height error, position error, the time
share of each state, and how many wind measurements the flight yields.
"""
import argparse, os, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import flightstates as xs

def w_colour(w):
    """Continuous colour from w: deep blue = rising, grey = calm, near-black
    = sinking. Clipped at +-2.5 m/s."""
    import numpy as np
    w = float(np.clip(w, -2.5, 2.5))
    grau = np.array([154, 150, 140.0])
    blau = np.array([18, 58, 140.0])
    dunkel = np.array([16, 20, 26.0])
    f = abs(w) / 2.5
    c = grau + (blau - grau) * f if w >= 0 else grau + (dunkel - grau) * f
    return tuple(c / 255.0)

COL_K = "#d68a00"          # circling   # sinking very strongly   # strongly sinking air   # straight, gliding

LW = {"K": 2.6, "G": 1.8}
BG = "#f2e4c8"
UTC_OFF = 2.0


def block_echt(df, t0, segs, ende):
    """The five rows from the full track."""
    n = len(df)
    std = (t0 + df["t"].to_numpy()) / 3600.0 + UTC_OFF
    alt = df["alt_g"].to_numpy()
    vario = df["vario"].to_numpy()
    dh, net = xs.turn_signal(df)
    _, _, wind = xs.wind_series(df, dh, net)
    circling = xs.circling_per_second(df, net, wind)
    kind = np.array([" "] * n, dtype="<U1")
    wsec = np.zeros(n)
    for s in segs:
        kind[s["a"]:s["b"] + 1] = s["kind"]
        if s["kind"] == "G":
            wsec[s["a"]:s["b"] + 1] = s["w"]
    kind[kind == " "] = segs[-1]["kind"]
    return dict(std=std, alt=alt, vario=vario, circling=circling, kind=kind,
                wsec=wsec, wind=wind, n=n)


def block_line(zl):
    """The same five rows, but from the text line alone."""
    _kennung, _glider, _pol, tag, segs, ende = xs.read_line(zl)
    t = np.array([s["t"] for s in segs] + [ende[0]], float)
    t = t + np.concatenate(([0], np.cumsum(np.diff(t) < 0) * 86400.0))
    h = np.array([s["h"] for s in segs] + [ende[1]], float)
    la = np.array([s["la"] for s in segs] + [ende[2]])
    lo = np.array([s["lo"] for s in segs] + [ende[3]])
    n = int(t[-1] - t[0]) + 1
    sek = np.arange(n) + t[0]
    alt = np.interp(sek, t, h)                       # eckig: Gerade je Segment
    vario = np.zeros(n)
    kind = np.array([" "] * n, dtype="<U1")
    wsec = np.zeros(n)
    wind = np.full(n, np.nan)
    wp_t, wp_v = [], []
    for i, s in enumerate(segs):
        a = int(t[i] - t[0]); b = int(t[i + 1] - t[0])
        kind[a:b] = s["kind"]
        vario[a:b] = (h[i + 1] - h[i]) / max(t[i + 1] - t[i], 1)
        if s["kind"] == "G":
            wsec[a:b] = s["w"]
        if s["kind"] in "Kk" and s.get("drift") is not None:
            wp_t.append(0.5 * (a + b)); wp_v.append(s["drift"])
    kind[kind == " "] = segs[-1]["kind"]
    if len(wp_t) >= 2:
        wind = np.interp(np.arange(n), wp_t, wp_v)
    elif wp_t:
        wind = np.full(n, wp_v[0])
    else:
        wind = np.zeros(n)
    circling = np.isin(kind, list("Kk"))
    std = sek / 3600.0 + UTC_OFF
    return dict(std=std, alt=alt, vario=vario, circling=circling, kind=kind,
                wsec=wsec, wind=wind, n=n), segs, ende, tag


def progress(std, kind, la=None, lo=None, segs=None, t0=None):
    """km/h straight-line per segment, as a staircase over time."""
    n = len(std)
    v = np.zeros(n)
    gr = np.concatenate(([0], np.flatnonzero(kind[1:] != kind[:-1]) + 1, [n]))
    for i in range(len(gr) - 1):
        a, b = gr[i], gr[i + 1] - 1
        dt = max(b - a, 1)
        dx = (lo[b] - lo[a]) * 111320.0 * np.cos(np.radians(la[a]))
        dy = (la[b] - la[a]) * 111132.0
        v[a:b + 1] = 3.6 * np.hypot(dx, dy) / dt
    return v


def draw(ax_alt, ax_var, ax_cir, ax_zus, ax_win, B, la, lo, title):
    std, kind, wsec = B["std"], B["kind"], B["wsec"]
    # 1 height
    ax_alt.plot(std, B["alt"], color="#a89f90", lw=0.6, zorder=2)
    wechsel = (kind[1:] != kind[:-1]) | (np.abs(wsec[1:] - wsec[:-1]) > 1e-9)
    gr = np.concatenate(([0], np.flatnonzero(wechsel) + 1, [len(kind)]))
    for i in range(len(gr) - 1):
        a, b = gr[i], gr[i + 1]
        if kind[a] == "K":
            c, dick = COL_K, LW["K"]
        else:
            c = w_colour(wsec[a]); dick = LW["G"] + 0.5 * min(abs(wsec[a]), 2.5)
        ax_alt.plot(std[a:b], B["alt"][a:b], color=c, lw=dick,
                    zorder=3, solid_capstyle="round")
    ax_alt.set_ylabel("altitude MSL [m]")
    ax_alt.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=10)
    # 2 vario
    v = B["vario"]
    ax_var.fill_between(std, 0, v, where=v >= 0, color=COL_K, alpha=.75, lw=0)
    ax_var.fill_between(std, 0, v, where=v < 0, color="#4a6a8a", alpha=.75, lw=0)
    ax_var.axhline(0, color="#333", lw=.6)
    ax_var.set_ylim(-4, 5); ax_var.set_ylabel("vario\n[m/s]")
    # 3 circling strip
    ax_cir.pcolormesh(std, [0, 1], B["circling"][np.newaxis, :len(std) - 1],
                      cmap=mcolors.ListedColormap(["#ffffff00", COL_K]),
                      vmin=0, vmax=1)
    ax_cir.set_yticks([]); ax_cir.set_ylabel("circling", rotation=0, ha="right",
                                             va="center")
    # 4 state strip
    farben = np.array([mcolors.to_rgb(COL_K) if k == "K" else w_colour(ww)
                       for k, ww in zip(kind, wsec)])
    ax_zus.imshow(farben[np.newaxis, :, :], aspect="auto",
                  extent=[std[0], std[-1], 0, 1])
    ax_zus.set_yticks([]); ax_zus.set_ylabel("state", rotation=0, ha="right",
                                             va="center")
    # 5 wind and progress
    vk = progress(std, kind, la, lo)
    ax_win.plot(std, vk, color="#333333", lw=0.9, label="progress per segment")
    ax_win.plot(std, B["wind"], color="#1f5fa8", lw=1.8,
                label="wind from circling drift")
    ax_win.set_ylim(0, 90); ax_win.set_ylabel("[km/h]")
    ax_win.legend(loc="upper left", fontsize=9, ncols=2, framealpha=.6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("igc")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.splitext(os.path.basename(a.igc))[0] + "_forensik.png"

    poldat = Path(__file__).resolve().parent / "polars.csv"
    polars = xs.load_polars(poldat) if poldat.exists() else None
    tag, glider, src, segs, ende, df = xs.segments(a.igc, polars)
    name = Path(a.igc).stem
    _, _, t0, _ = xs.read_igc(a.igc)
    zl = xs.line("P", glider, src, tag, segs, ende)

    Be = block_echt(df, t0, segs, ende)
    Br, segs_r, ende_r, _ = block_line(zl)

    la_e, lo_e = df["lat"].to_numpy(), df["lon"].to_numpy()
    # reconstructed positions: straight lines between the segment starts
    t = np.array([s["t"] for s in segs_r] + [ende_r[0]], float)
    t = t + np.concatenate(([0], np.cumsum(np.diff(t) < 0) * 86400.0))
    la_r = np.interp(np.arange(Br["n"]) + t[0], t,
                     [s["la"] for s in segs_r] + [ende_r[2]])
    lo_r = np.interp(np.arange(Br["n"]) + t[0], t,
                     [s["lo"] for s in segs_r] + [ende_r[3]])

    W, H, DPI = 9000, 2600, 100
    fig, ax = plt.subplots(11, 1, figsize=(W / DPI, H / DPI), dpi=DPI,
                           sharex=True, gridspec_kw=dict(
                               height_ratios=[5, 1.7, .6, .6, 1.5, 1.0,
                                              5, 1.7, .6, .6, 1.5],
                               hspace=.07))
    fig.patch.set_facecolor("#eef0ea")
    ax[5].set_axis_off()
    ax = np.concatenate([ax[:5], ax[6:]])
    for x in ax:
        x.set_facecolor(BG); x.grid(True, axis="x", color="#c9b697", lw=.6)
        x.margins(x=.001)

    kb = len(zl.encode())
    draw(*ax[0:5], Be, la_e, lo_e,
            f"TOP — from the real IGC file, {Be['n']} points at 1 Hz   ·   {name}   ·   "
            f"{len(segs)} segments   ·   colours: "
            + "colour = w, the movement of the air; orange = circling")
    draw(*ax[5:10], Br, la_r, lo_r,
            f"BOTTOM — from the text line alone, {len(segs_r)} "
            f"segment starts, {kb} characters   ·   per segment only kind, second, "
            f"height, position; circling segments also carry turn count and drift")

    ticks = np.arange(np.ceil(Be["std"][0] * 4) / 4, Be["std"][-1], .25)
    for x in ax:
        x.set_xticks(ticks)
    ax[-1].set_xticklabels([f"{int(t):02d}:{int(round((t % 1) * 60)):02d}"
                            if abs((t * 2) % 1) < 1e-6 else "" for t in ticks],
                           fontsize=9)
    ax[-1].set_xlabel("local time — one grid line per 15 min")
    fig.savefig(out, facecolor=fig.get_facecolor())

    # ---- numeric comparison: original against reconstruction -------------
    m = min(Be["n"], Br["n"])
    dh_ = Br["alt"][:m] - Be["alt"][:m]
    dx = (lo_r[:m] - lo_e[:m]) * 111320.0 * np.cos(np.radians(la_e[:m]))
    dy = (la_r[:m] - la_e[:m]) * 111132.0
    do = np.hypot(dx, dy)
    print(f"  height   : mean error {np.abs(dh_).mean():5.1f} m, "
          f"90% under {np.percentile(np.abs(dh_),90):5.1f} m, "
          f"largest {np.abs(dh_).max():5.0f} m")
    print(f"  position : mean error {do.mean():5.0f} m, "
          f"90% under {np.percentile(do,90):5.0f} m, "
          f"largest {do.max():5.0f} m")
    print("  time share per state (real / reconstructed):")
    NAME = {"K": "circling", "G": "straight"}
    for k in "KG":
        pe = 100 * (Be["kind"][:m] == k).mean(); pr = 100 * (Br["kind"][:m] == k).mean()
        print(f"    {NAME[k]:22s} {pe:5.1f} %  /  {pr:5.1f} %")
    we = Be["wsec"][:m][Be["kind"][:m] == "G"]
    wr = Br["wsec"][:m][Br["kind"][:m] == "G"]
    print(f"    w over the straight seconds: real mean {we.mean():+.2f}, "
          f"reconstructed {wr.mean():+.2f}")
    w = [s2["drift"] for s2 in segs if s2.get("drift") is not None]
    print(f"  wind     : {len(w)} measurements from "
          f"{sum(1 for s2 in segs if s2['kind'] in 'Kk')} circling segments, "
          f"mean {np.mean(w):.1f} km/h, largest {np.max(w):.1f} km/h")
    print(f"{out}  {W}x{H}  line {kb} characters, {len(segs)} segments")
    return zl


if __name__ == "__main__":
    main()
