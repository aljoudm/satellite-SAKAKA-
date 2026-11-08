# GEO Ku-band link for rural broadband in Al-Jouf (Sakaka) — Section 6: Modulation and Multiple Access

A GEO satellite at 45°E serves rural broadband users in Al-Jouf, Saudi Arabia (reference site: Sakaka). The forward link uses DVB-S2 8PSK with code rate 3/4 (spectral efficiency 2.228 bit/symbol) at 30 Msym/s with a root-raised-cosine filter (roll-off 0.20), which fits the 36 MHz transponder. The overall Es/N0 is 17.35 dB in clear sky, 16.44 dB at 99.9 % availability and 14.72 dB at 99.99 %, which leaves 7.0 dB of margin at 99.9 % after a 1.5 dB implementation margin. The return link uses MF-TDMA with 14 carriers of 2 Msym/s (QPSK 3/4, 2.97 Mbps per carrier, 41.6 Mbps total), shared by demand assignment.

## How to run

Requires MATLAB (tested on R2026a) with Communications Toolbox and Signal Processing Toolbox. From MATLAB:

```matlab
run('section6_modulation_access.m')
```

or from a terminal: `matlab -batch "run('section6_modulation_access.m')"` (about 2 minutes). Figures go to `figures/`, numbers to `results/summary.txt` and `results/results.mat`. Set `runLDPC = true` at the top of the script to add the optional DVB-S2 LDPC frame-error-rate check (off by default; not exercised in the delivered run).

## Figures

| File | Content |
|---|---|
| `figures/Fig6_1_constellation.png` | Received Gray-coded QPSK and 8PSK constellations at Es/N0 = 16.44 dB |
| `figures/Fig6_2_BER.png` | Uncoded BER vs Eb/N0, theory and Monte-Carlo, with the 1e-6 line |
| `figures/Fig6_3_shannon.png` | Spectral efficiency vs Eb/N0: Shannon bound and the DVB-S2 MODCODs |
| `figures/Fig6_4_margin.png` | Required Es/N0 per MODCOD and margin at 99.9 % / 99.99 % availability |
| `figures/Fig6_5_spectrum.png` | Transmitted 8PSK spectrum (RRC 0.20, 30 Msym/s) with transponder edges |
| `figures/Fig6_6_mftdma.png` | MF-TDMA example frame and histogram of active homes per frame |

## Summary results

| MODCOD | η (bit/sym) | Ideal Es/N0 (dB) | Required (dB) | Gap to Shannon (dB) | Margin 99.9 % (dB) | Margin 99.99 % (dB) |
|---|---|---|---|---|---|---|
| QPSK 1/2 | 0.989 | 1.00 | 2.50 | 1.07 | 13.94 | 12.22 |
| QPSK 3/4 | 1.487 | 4.03 | 5.53 | 1.47 | 10.91 | 9.19 |
| 8PSK 2/3 | 1.981 | 6.62 | 8.12 | 1.93 | 8.32 | 6.60 |
| **8PSK 3/4** | **2.228** | **7.91** | **9.41** | **2.25** | **7.03** | **5.31** |
| 8PSK 5/6 | 2.479 | 9.35 | 10.85 | 2.75 | 5.59 | 3.87 |

- Uncoded BER = 1e-6 at Eb/N0 = 10.53 dB (QPSK) and 13.95 dB (8PSK).
- Shannon minimum Eb/N0 at η = 2.228: 2.19 dB.
- 99 % occupied bandwidth 32.2 MHz (nominal Rs(1+α) = 36 MHz); 99.98 % of the power lies inside the ±18 MHz transponder.
- MF-TDMA (300 homes, p = 0.05, 2 Mbps each, 10,000 frames): mean demand 29.9 Mbps (15.0 active homes), mean utilisation 71.9 % offered / 71.3 % carried, demand exceeds the 41.6 Mbps capacity in 7.7 % of frames, 99th-percentile demand 48 Mbps (24 homes). Fixed FDMA would need 300 × 2 = 600 Mbps, 14.4 times the shared capacity.
