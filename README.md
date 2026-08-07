# Multi-target tracking of single targets in broadband split-beam echosounder data

MATLAB code and example data for detecting and tracking individual acoustic
targets in broadband split-beam recordings, using a converted measurement
Kalman filter (CMKF).

Each ping is pulse compressed and searched for echo peaks; the split-aperture
correlator gives the angle of arrival of every accepted peak. The resulting
range/angle measurements are converted to Cartesian coordinates before
filtering, so the filter itself stays linear and only the measurement
covariance carries the nonlinearity of the conversion. Its Jacobian is
evaluated at the predicted position, which keeps the converted covariance
independent of the current measurement noise.

The repository reproduces the detection and tracking figure in the paper from
a five-minute recording made by an upward-looking, moored WBAT at 590 m depth,
and is set up to be pointed at other recordings.

## Quick start

```matlab
run_example
```

The script adds `src/` and `external/` to the path, processes
`data/example_pings_20211105.mat`, prints a summary, and draws the five-panel
overview figure. It takes well under a minute.

Expected output:

```
150 pings, 05-Nov-2021 10:30:14 to 05-Nov-2021 10:35:18, sound speed 1487.8 m/s
Detections: 6270
Tracking completed in 1.2 s
Confirmed tracks: 205 (63 still active at the last ping)
Tracks with more than 10 hits: 145
Fitted speed [m/s] - min 0.016 | median 0.071 | max 0.204
Fitted R^2 median [x y z] - [0.752 0.804 0.596]
```

## Requirements

- MATLAB **R2023b or later** (`orderedcolors`, `xsecondarylabel`)
- **Signal Processing Toolbox** (`resample`, `findpeaks`, `tukeywin`, `db`)

No other toolboxes are needed. `gsw_sound_speed` and `viridis` are bundled in
`external/`.

## Layout

```
run_example.m                 Driver: reproduces the figure from the paper
src/
  defaultConfig.m             Every tunable parameter, in one struct
  loadPings.m                 WBAT/EK80 export -> compact ping struct
  soundSpeedFromCTD.m         Sound speed at the ping times, from moored CTD
  detectTargets.m             Peak detection, angle of arrival, target strength
  trackTargetsCMKF.m          Multi-target CMKF tracker
  cartesianToSplitBeam.m      [x y z] -> [r alpha beta]
  convertedMeasurementCovariance.m
                              Measurement covariance through the conversion
  constantVelocityModel.m     State transition F and process noise Q
  fitTrackVelocity.m          Least-squares velocity per track
  plotTrackingOverview.m      The five-panel figure
data/                         Example recording, transmit signal, CTD casts
external/                     Third-party code, see external/README.md
tools/
  makeExampleDataset.m        How the shipped data/ files were reduced
```

## Data flow

```
loadPings ─┐
           ├─> detectTargets ──> detections (table)  ──> trackTargetsCMKF ──> tracks
soundSpeed ┘                     env (geometry)      ──┘        │
                                                                v
                                                        fitTrackVelocity
                                                                │
                                                                v
                                                       plotTrackingOverview
```

`detections` is a table with one row per accepted echo peak; `env` carries the
range axis, ping times and sampling geometry that the tracker and the figure
both need; `tracks` is a struct array of confirmed tracks, each holding its
state, covariance, position log and the detection rows assigned to it. The
field-by-field meaning of each is in the header of the function that produces
it (`help detectTargets`, `help trackTargetsCMKF`).

## Running on your own data

1. **Reading your files.** `loadPings` expects the WBAT export layout described
   in its header. For any other format, write a replacement that returns the
   same fields — nothing downstream reads the raw file.

2. **Sound speed.** `soundSpeedFromCTD` interpolates the moored CTD records in
   `data/`. With no CTD record, pass a constant to `detectTargets` instead.

3. **Parameters.** Everything tunable lives in `defaultConfig`, split into
   `cfg.detection`, `cfg.tracking` and `cfg.analysis`. The values that most
   often need changing for a different deployment are:

   | Parameter | Meaning |
   | --- | --- |
   | `detection.maxAnalysisRange` | how far out to process [m] |
   | `detection.minDetectionRange` | near-field blanking [m] |
   | `detection.powerThreshold` | detection floor [linear power] |
   | `detection.peakWidthLimits` | expected echo length, set by the pulse [m] |
   | `detection.transducerRadius` | effective aperture radius of your transducer [m] |
   | `tracking.sigmaR`, `tracking.kappaTheta` | measurement noise, range and angular |
   | `tracking.sigmaAccel` | how hard targets are expected to manoeuvre |
   | `analysis.instrumentDepth` | transducer depth, for the depth axes [m] |

   Detection windows quoted in metres are rescaled by the in situ sound speed
   at run time, so they keep a constant extent as sound speed changes.

4. **Geometry.** The angle scaling assumes a three-sector split-beam
   transducer whose sectors are 120-degree pie slices, and the beam-pattern
   compensation assumes a circular piston. Both are set by
   `detection.transducerRadius`.

## Example data

`data/` holds reduced copies of the original recordings, small enough to
distribute with the code:

| File | Contents | Reduced from |
| --- | --- | --- |
| `example_pings_20211105.mat` | 150 pings, three sectors, pulse compressed | `AR59_2-Phase4-D20211105-T103012-0.mat` |
| `transmit_signal.mat` | transmit signal `ytx` and its rate `fs_dec` | `WBAT_tx.mat` |
| `ctd_2021.mat`, `ctd_2022.mat` | moored CTD time series | `WBAT_CTD_2021.mat`, `WBAT_CTD_2022.mat` |

Two reductions were applied, both to keep the repository small and neither
affecting the results:

- the raw, uncompressed samples (`complexsamples`) were dropped, since only the
  pulse-compressed samples are used, and the pulse-compressed samples were
  truncated to the top 64 m — the analysis reaches 45 m;
- the CTD records keep only the four fields used to compute sound speed.

`tools/makeExampleDataset.m` performs both steps and documents the provenance.
Setting `cfg.detection.maxAnalysisRange` beyond 64 m raises an error against
the shipped file; use an untruncated export to go deeper.

## Citation

See `CITATION.cff`.

## Licence

MIT, see `LICENSE`. Third-party code in `external/` keeps its own licence; see
`external/README.md`.
