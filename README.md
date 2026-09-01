# flightstates — flight track in, one line of text out

Reads IGC tracks. Writes one text line per flight: the flight as a sequence
of segments. Two kinds, everything else is a number:

| Code | Segment | carries |
|---|---|---|
| `K` | circling | turn count, drift = **the wind** |
| `G` | straight | `w` = **vertical movement of the air** (m/s, own sink removed), airspeed, path |

The unit is the second, not the thermal. So the air that carries a pilot
*without* circling — the air that decides a cross-country flight — is in the
data, not lost between climbs. No classes: maps colour continuously from `w`.

`flightstates.py` runs on its own; needs `numpy`, `pandas`, `scipy`. No maps,
no weather, no network.

## The runs

Removing the glider's own sink from the vario needs a sink polar per glider
type. It is measured first and inspected before the segmenter runs. On a
big or foreign archive, measure a sample first and look at it:

```bash
pip install numpy pandas scipy

python3 polarmaker.py "flights/**/*.IGC" --part 1/10 --out probe.npz   # 1  sample
python3 polarmaker.py --join probe.npz --out probe.csv --min 5         #    look at probe.csv + probe_stat.csv
python3 polarmaker.py "flights/**/*.IGC" --out polars.csv --min 50     # 2  the polar table
python3 flightstates.py --delta "flights/**/*.IGC" > states.txt        # 3  the segmented flights
```

On a small archive of your own, step 1 can be skipped. What the pieces mean:

* `"flights/**/*.IGC"` — replace `flights` with the folder holding your IGC
  files; `**` searches all subfolders. The quotes are required (the script
  expands the pattern itself; case of `.igc` doesn't matter). Several
  patterns may be given. Duplicate files (same name and size) are
  processed once.
* `--min 50` — a glider type gets its own row only with at least 50 flights.
* `--delta` — compact numbers (recommended); without it, plain numbers.
* `> states.txt` — the segmenter prints to standard output; `>` puts it in
  a file.
* `--polars FILE` (step 3) — which table to use; default: `polars.csv` next
  to the script. `--id TEXT` — what to write in the first field (default `P`).

`polars.csv` is one small table: per glider type, own sink at 25–65 km/h,
measured as the mode of the vario per band (the most common air on a glide
is near-still air). Values that measured air instead of glider — slow
flight spent in ridge lift does this — are discarded by physics and
stability checks (thresholds and their justification: header of
`polarmaker.py`). **No positions, no times, no pilots** — every line
readable by eye. Glider type comes from `HFGTY`; names are normalised but
digits are never touched. Per single flight the sink estimate scatters
0.21 m/s; per type at 50+ flights, ±0.03. Flights without a table row use
the `_general` row; the line marks it (`pol` field: `g` own row, `a`
general).

Next to the table, `polars_stat.csv` records **every band of every glider**
— kept or discarded, with the reason and a 90 % confidence interval from a
bootstrap over flights. Which glider types to trust (or blacklist) can be
decided later from this file and the `.npz` sample dumps alone, without
touching the IGC files again.

Both runs parallelise with `--part k/n`. Polar parts write `.npz`;
`--join "teil*.npz"` builds the table, identical to a single pass.
Segmenter parts are plain text; `cat` them. Memory stays small: per flight
only a 4 KB histogram is kept, never the raw seconds.

## Running it for someone else

If you hold an archive and someone asks you to compress it for them, the
steps above are all there is to do, and **two files** come out:

| | |
|---|---|
| `polars.csv` | one small table, a few hundred bytes — one line per glider type |
| `states.txt` | one line per flight; with `--delta` and `xz -6` about 364 bytes per flying hour |

Send those two. Nothing else is needed, and nothing else should be sent: the
IGC files stay where they are. Every line is readable by eye before it leaves
the house — no pilot names, and in `polars.csv` no positions or times at all.
The glider type is the only thing carried over from the IGC header; it is
what the sink polar is fitted per, and without it every flight falls back to
the `_general` row.

`polars.csv` belongs with the shipment: the header line of `states.txt`
records its SHA1, so the recipient can check which table was used.

## The output line

```
# flightstates 2.0 polars=polars.csv sha1=a48b1d37a814      <- file header, once
id;glider;pol;yyyymmdd;SEG|SEG|...;END

G (straight):   G,t,h,lat,lon,w,v,z
K (circling):   K,t,h,lat,lon,turns,drift_kmh,drift_from_deg
END:              t,h,lat,lon
```

`t,h,lat,lon` are second of day (UTC), altitude MSL and position at the
segment's **beginning**; the end of one segment is the beginning of the next.
`w` = vertical air movement in m/s, own sink removed. `v` = true airspeed
km/h (mean over the seconds). `z` = flown path over chord, in percent.
With `v` and the height difference, a sink polar can be re-derived from the
lines alone. `pol` records where the own-sink curve came from: `g` glider's
row, `a` `_general`, `n` none (constant −1.11).

**The first field is yours.** Default is the letter `P`; `--id TEXT` sets it.
An anonymised pilot key would make technique-over-the-years studies possible —
an offer, not a request; nothing needs the field.

`--delta` writes every number as difference to the previous — same content,
short numbers; `read_line()` / `read_line_delta()` in `flightstates.py` turn
a line back into segments.

```
P;NIVIUK Artik 6;g;20230522;k27039,1920,4637603,802630,0.6,,|R45,-3,-160,42,-0.2|…
```

## How it decides

Per second: circling if the summed heading change over a centred 20 s window
exceeds 60° (checked against hand-marked flights: 94.7 % agreement). Wind per
thermal from the circling drift (≥1.5 turns, ≥30 s, endpoints trimmed,
≤40 km/h), interpolated between thermals. Straight stretches are cut on the
**air-corrected height** (height minus accumulated own sink at the flown
airspeed) and the ground plan — joint Douglas-Peucker, 30 m height, 60 m
position. Boundaries are pure geometry; no label enters the cutting.

All thresholds are constants at the top of `flightstates.py`:

| | | |
|---|---|---|
| sampling | 1 Hz | turn window/threshold | 20 s / 60° |
| height smoothing | Savitzky-Golay 15 s | wind sample | ≥1.5 turns, ≥30 s, ≤40 km/h |
| tolerances | 30 m height, 60 m position | wind smoothing | 600 s |
| minimum run | 20 s | recorder gaps >60 s | longest stretch kept |
| coordinates | 5 decimals (~1 m) | altitude spikes >500 m | dropped |

## Check

```bash
python3 verify.py flight.IGC --out check.png
python3 chart.py  flight.IGC --out plate.png
```

`verify.py` draws the flight twice — from the IGC, and from the text line
alone — and prints the errors. Typical: height mean 13 m (90 % under 31 m),
position mean 39 m (90 % under 88 m), every `w` and all time shares identical.

`chart.py` draws one flight as a plate like the picture below. Options:
`--raw` adds the untouched trace on top, `--language en`, `--title TEXT`.
Both scripts expect `polars.csv` and `flightstates.py` beside them.

![one flying day: IGC above, the line below](example_day.png)

## Size and time

Measured on 774 alpine flights (2 883 h): **9.2 bytes per segment, 364 bytes
per flying hour** (`--delta`, `xz -6`); both runs together ~30 ms per flying
hour per core. A world archive of ~950 000 flights: **1–2 GB compressed,
about one core-day** — an afternoon on eight cores. Lines are independent;
process the file line by line.

```bash
for k in 1 2 3 4; do
  python3 flightstates.py --delta --part=$k/4 "flights/**/*.IGC" > part$k.txt &
done
wait; cat part*.txt | xz -6 > states.txt.xz
```

## Files

| | |
|---|---|
| `flightstates.py` | run 2 — the segmenter, self-contained |
| `polarmaker.py` | run 1 — sink polar per glider type |
| `polars.csv` | the shipped table (774 flights); replace with your own |
| `verify.py`, `chart.py` | check picture, day plate |
| `example_line.txt`, `example_line_delta.txt` | one flight, both forms |

Notes: barometric height when plausible, else GPS; date from `HFDTE`, times
UTC; flights without valid B records are skipped and reported on stderr. The
wind method is taken unchanged from a tested implementation. The polar is
measured in near-still air; it contains harness and trim, and is not a
manufacturer polar.

Jörg Korner, 2026. Free to use.
