#!/usr/bin/env python3
"""chart.py — one plate: a flying day in seven states.

  python3 chart.py flight.IGC --out plate.png --language en
  python3 chart.py flight.IGC --out plate.png --raw

--raw puts the untouched trace on top, exactly as it comes out of the IGC
file, and the same trace below it with every second sorted into a state. That
is the whole difference this package makes, in one picture.
"""
import argparse, os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
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

COL_K = "#d68a00"          # circling   # sinking very strongly   # strongly sinking air
LW = {"K": 2.6, "G": 1.8}
BG = "#f2e4c8"
W = {"de": dict(K="Kreisen", hoehe="Höhe MSL [m]", vario="Vario\n[m/s]",
                zust="Zustand", wind="[km/h]", zeit="Ortszeit",
                wl="Wind aus der Kreisdrift", vk="Vorankommen je Segment",
                roh="Die IGC-Datei, so wie sie ankommt",
                unser="Dieselbe Spur, jede Sekunde einem Zustand zugeordnet"),
     "en": dict(K="circling", hoehe="altitude MSL [m]", vario="vario\n[m/s]",
                zust="state", wind="[km/h]", zeit="local time",
                wl="wind from circling drift", vk="progress per segment",
                roh="The IGC file, as it arrives",
                unser="The same track, every second sorted into a state")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("igc"); ap.add_argument("--out", default="tafel.png")
    ap.add_argument("--language", default="de", choices=("de", "en"))
    ap.add_argument("--width", type=int, default=4500)
    ap.add_argument("--raw", action="store_true",
                    help="add the untouched trace as a panel on top")
    ap.add_argument("--title", default="")
    a = ap.parse_args(); w = W[a.language]

    poldat = Path(__file__).resolve().parent / "polars.csv"
    polars = xs.load_polars(poldat) if poldat.exists() else None
    tag, glider, src, segs, ende, df = xs.segments(a.igc, polars)
    _, _, t0, _ = xs.read_igc(a.igc)
    n = len(df)
    std = (t0 + df["t"].to_numpy()) / 3600.0 + 2.0
    alt = df["alt_g"].to_numpy(); vario = df["vario"].to_numpy()
    la = df["lat"].to_numpy(); lo = df["lon"].to_numpy()
    kind = np.array([" "] * n, dtype="<U1")
    wsec = np.zeros(n)
    for s in segs:
        kind[s["a"]:s["b"] + 1] = s["kind"]
        if s["kind"] == "G":
            wsec[s["a"]:s["b"] + 1] = s["w"]
    kind[kind == " "] = segs[-1]["kind"]
    wind = np.full(n, np.nan); pt, pv = [], []
    for s in segs:
        if s["kind"] in "Kk" and s.get("drift") is not None:
            pt.append(0.5 * (s["a"] + s["b"])); pv.append(s["drift"])
    wind = np.interp(np.arange(n), pt, pv) if len(pt) >= 2 else np.zeros(n)
    vk = np.zeros(n)
    for s in segs:
        d = max(s["b"] - s["a"], 1)
        vk[s["a"]:s["b"] + 1] = 3.6 * np.hypot(
            (lo[s["b"]] - lo[s["a"]]) * 111320 * np.cos(np.radians(la[s["a"]])),
            (la[s["b"]] - la[s["a"]]) * 111132) / d

    DPI = 100
    hr = ([4.2, 5, 1.6, .7, 1.2, 1.2] if a.raw else [5, 1.6, .7, 1.2, 1.2])
    H = int(a.width * (0.225 if a.raw else 0.175))
    fig, ax = plt.subplots(len(hr), 1, figsize=(a.width / DPI, H / DPI), dpi=DPI,
                           sharex=True, gridspec_kw=dict(height_ratios=hr, hspace=.16 if a.raw else .08))
    fig.patch.set_facecolor("#eef0ea")
    for x in ax:
        x.set_facecolor(BG); x.grid(True, axis="x", color="#c9b697", lw=.6); x.margins(x=.002)
    if a.raw:
        roh = ax[0]
        roh.plot(std, alt, color="#3a3833", lw=1.0)
        roh.set_ylabel(w["hoehe"])
        roh.set_title(w["roh"], loc="left", fontsize=12, fontweight="bold", pad=8)
        ax = ax[1:]
        ax[0].set_title(w["unser"], loc="left", fontsize=12, fontweight="bold", pad=8)
    ax[0].plot(std, alt, color="#a89f90", lw=.6, zorder=2)
    for s_ in segs:
        p, q = s_["a"], min(s_["b"] + 1, n - 1) + 1
        if s_["kind"] == "K":
            farbe, dick = COL_K, LW["K"]
        else:
            farbe, dick = w_colour(s_["w"]), LW["G"] + 0.5 * min(abs(s_["w"]), 2.5)
        ax[0].plot(std[p:q], alt[p:q], color=farbe, lw=dick,
                   zorder=3, solid_capstyle="round")
    ax[0].set_ylabel(w["hoehe"])
    from matplotlib.lines import Line2D
    hand = [Line2D([], [], color=COL_K, lw=3, label=w["K"])]
    for wv in (2.0, 1.0, 0.0, -1.0, -2.0):
        hand.append(Line2D([], [], color=w_colour(wv), lw=3,
                           label=f"w {wv:+.0f} m/s"))
    ax[0].legend(handles=hand, loc="upper left", ncols=6, fontsize=11,
                 framealpha=.7)
    if a.title and not a.raw:
        ax[0].set_title(a.title, loc="left", fontsize=12, fontweight="bold", pad=8)
    elif a.title and a.raw:
        fig.suptitle(a.title, x=0.008, ha="left", fontsize=13, fontweight="bold")
    ax[1].fill_between(std, 0, vario, where=vario >= 0, color=COL_K, alpha=.75, lw=0)
    ax[1].fill_between(std, 0, vario, where=vario < 0, color="#4a6a8a", alpha=.75, lw=0)
    ax[1].axhline(0, color="#333", lw=.6); ax[1].set_ylim(-4, 5); ax[1].set_ylabel(w["vario"])
    farben = np.array([mcolors.to_rgb(COL_K) if k == "K" else w_colour(ww)
                       for k, ww in zip(kind, wsec)])
    ax[2].imshow(farben[np.newaxis, :, :], aspect="auto",
                 extent=[std[0], std[-1], 0, 1])
    ax[2].set_yticks([]); ax[2].set_ylabel(w["zust"], rotation=0, ha="right", va="center")
    ax[3].plot(std, vk, color="#333", lw=.9, label=w["vk"])
    ax[3].set_ylim(0, 90); ax[3].set_ylabel(w["wind"])
    ax[3].legend(loc="upper left", fontsize=10, framealpha=.6)
    # the wind on its own scale, one dot per sample: a long climb carries
    # several (one per piece of ~4 turns) — the profile through the climb
    ax[4].plot(std, wind, color="#1f5fa8", lw=1.2, label=w["wl"])
    if len(pt) >= 2:
        ax[4].plot(np.interp(pt, np.arange(n), std), pv, "o", ms=3.6,
                   color="#1f5fa8", mec="white", mew=.6)
    ax[4].set_ylim(0, max(10.0, float(np.nanmax(pv)) * 1.25 if len(pv) else 10.0))
    ax[4].set_ylabel(w["wind"])
    ax[4].legend(loc="upper left", fontsize=10, framealpha=.6)
    ticks = np.arange(np.ceil(std[0] * 4) / 4, std[-1], .25)
    for x in ax:
        x.set_xticks(ticks)
    ax[4].set_xticklabels([f"{int(t):02d}:{int(round((t % 1) * 60)):02d}"
                           if abs((t * 2) % 1) < 1e-6 else "" for t in ticks], fontsize=10)
    ax[4].set_xlabel(w["zeit"])
    fig.savefig(a.out, facecolor=fig.get_facecolor())
    print(a.out, a.width, H, len(segs), "segments")


if __name__ == "__main__":
    main()
