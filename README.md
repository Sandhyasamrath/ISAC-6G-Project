# Optimal Joint Radar-Communication Waveform Design for 6G ISAC

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![6G](https://img.shields.io/badge/6G-ISAC-red.svg)](#)
[![mmWave](https://img.shields.io/badge/mmWave-28%20GHz-green.svg)](#)

> A complete simulation framework for **Integrated Sensing and Communication (ISAC)** in 6G networks, demonstrating the trade-off between communication throughput and radar sensing accuracy in an OFDM-based mmWave MIMO system.

---

## 📋 Abstract

This project presents an end-to-end MATLAB simulation framework for an OFDM-based ISAC system operating in the **28 GHz mmWave band** with **100 MHz bandwidth** and **8×4 MIMO** configuration. The framework implements joint waveform design where a single OFDM signal simultaneously serves communication users and detects radar targets, characterizing the fundamental performance trade-off between sensing accuracy and communication throughput through a **Pareto optimization framework**.

---

## 🎯 Problem Statement

The central problem is the **fundamental performance trade-off between communication and sensing** in an ISAC system. Both functions share:
- Same OFDM waveform
- Same hardware (antennas, RF chains)
- Same time-frequency resources
- Same transmit power

Any design choice that improves one function tends to degrade the other. This project characterizes this trade-off using two key parameters:

| Parameter | Symbol | Effect |
|-----------|--------|--------|
| **Pilot subcarrier ratio** | α | More pilots → better sensing, lower throughput |
| **ISAC weight** | λ | Higher λ → communication priority; lower λ → sensing priority |

---

## 🚀 Key Results

All results obtained at **15 dB operating SNR**:

### 📡 Radar Performance
| Metric | Value |
|--------|-------|
| Radar SINR | **18.5 dB** |
| Detected targets | **3 / 3** ✅ |
| Range RMSE | **1.16 m** |
| Range resolution | **1.92 m** |

### 📶 Communication Performance
| Metric | Value |
|--------|-------|
| Spectral efficiency | **5.03 bps/Hz** |
| Bit error rate (16-QAM) | **~10⁻³** |
| Total bits transmitted | 74,400 |

### 📐 Theoretical Bounds (CRB)
| Metric | Value |
|--------|-------|
| CRB Range | **0.22 m** |
| CRB Velocity | **0.86 m/s** |
| Effective bandwidth | 22.9 MHz |

### ⚖️ Pareto Optimal Operating Point
| Parameter | Optimal Value |
|-----------|---------------|
| **α\*** (pilot ratio) | **0.20** |
| **λ\*** (ISAC weight) | **0.50** |
| Communication SE | 75% of pure-comm peak |
| Sensing Quality | 70% of pure-sensing peak |

---

## 🏗️ System Architecture

The framework consists of **9 modular sections**:

```
┌─────────────────────────────────────────────────────────┐
│  S1 → S2 → S3 → S4 → S5 → S6 → S7 → S8 → S9             │
└─────────────────────────────────────────────────────────┘
   │     │     │     │     │     │     │     │     │
   │     │     │     │     │     │     │     │     └─ Final Visualization
   │     │     │     │     │     │     │     └─────── Pareto Surface U(λ,α)
   │     │     │     │     │     │     └───────────── CRB & Detection Probability
   │     │     │     │     │     └─────────────────── Communication Receiver
   │     │     │     │     └───────────────────────── Radar Receiver
   │     │     │     └─────────────────────────────── Dual-Path Channel
   │     │     └───────────────────────────────────── OFDM-ISAC Transmitter
   │     └─────────────────────────────────────────── MIMO Beamforming
   └───────────────────────────────────────────────── System Parameters
```

---

## ⚙️ System Parameters

| Parameter | Value |
|-----------|-------|
| Carrier frequency | 28 GHz (mmWave) |
| Bandwidth | 100 MHz |
| Total subcarriers | 256 (200 active + 56 guard) |
| OFDM symbols per CPI | 128 (4 pilot + 124 data) |
| Cyclic prefix | 32 samples |
| Modulation | 16-QAM |
| Transmit antennas | 8 (ULA) |
| Receive antennas | 4 (ULA) |
| Antenna spacing | λ/2 |
| Operating SNR | 15 dB |

### 🎯 Radar Targets Configuration
| Target | Range | Velocity | RCS | Angle |
|--------|-------|----------|-----|-------|
| T1 | 80 m | +15 m/s | 1.0 | +10° |
| T2 | 150 m | −20 m/s | 0.6 | −20° |
| T3 | 280 m | +30 m/s | 0.8 | +35° |

---

## 🚀 How to Run

### Prerequisites
- MATLAB R2021a or later
- Signal Processing Toolbox
- Communications Toolbox (recommended)

### Steps

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Sandhyasamrath/ISAC-6G-Project.git
   cd ISAC-6G-Project/code/matlab
   ```

2. **Open MATLAB** and navigate to the `code/matlab/` directory

3. **Run the master script:**
   ```matlab
   run_all
   ```

   This will execute all 9 modules sequentially and produce all plots.

4. **Or run individual modules** in order:
   ```matlab
   ISAC_1   % System parameters
   ISAC_2   % MIMO beamforming
   ISAC_3   % OFDM-ISAC transmitter
   ISAC_4   % Dual-path channel
   ISAC_5   % Radar receiver
   ISAC_6   % Communication receiver
   ISAC_7   % CRB analysis
   ISAC_8   % Pareto surface
   ISAC_9   % Final visualization
   ```

   ⚠️ Each script saves a workspace `.mat` file consumed by the next, so run in order.

---

## 📁 Project Structure

```
ISAC-6G-Project/
├── README.md
├── code/
│   └── matlab/
│       ├── ISAC_1.m  →  System parameters
│       ├── ISAC_2.m  →  MIMO beamforming
│       ├── ISAC_3.m  →  OFDM-ISAC transmitter
│       ├── ISAC_4.m  →  Dual-path channel
│       ├── ISAC_5.m  →  Radar receiver (Range-Doppler + CFAR)
│       ├── ISAC_6.m  →  Communication receiver (LS + ZF)
│       ├── ISAC_7.m  →  CRB analysis
│       ├── ISAC_8.m  →  Pareto surface optimization
│       ├── ISAC_9.m  →  Final visualization
│       └── run_all.m →  Master script
├── docs/
│   ├── ISAC_Project_Report.docx       →  Full project report
│   └── ISAC_Project_Presentation.pptx →  Project presentation
└── references/
    ├── Optimal-Joint-ISAC-Waveform-Design....pdf
    └── Optimized_Waveform_Design_for_OFDM-based_ISAC.pdf
```

---

## 📚 References

1. **C. Sturm and W. Wiesbeck**, "Waveform design and signal processing aspects for fusion of wireless communications and radar sensing," *Proc. IEEE*, vol. 99, no. 7, pp. 1236–1259, Jul. 2011.

2. **J. A. Zhang et al.**, "An overview of signal processing techniques for joint communication and radar sensing," *IEEE J. Sel. Topics Signal Process.*, vol. 15, no. 6, pp. 1295–1315, Nov. 2021.

3. **F. Liu et al.**, "Cramér-Rao Bound optimization for joint radar-communication beamforming," *IEEE Trans. Signal Process.*, vol. 70, pp. 240–253, Jan. 2022.

4. **J. Singh et al.**, "Pareto optimal hybrid beamforming for short-packet mmWave ISAC," *arXiv:2406.01945*, Jun. 2024.

5. **S. Mura et al.**, "Optimized waveform design for OFDM-based ISAC under limited resource occupancy," *arXiv:2406.19036*, Jun. 2024.

---

## 🎓 Project Information

- **Title:** Optimal Joint Radar-Communication Waveform Design and Performance Trade-Offs in 6G ISAC
- **Domain:** Wireless Communication, Signal Processing, 6G Networks
- **Year:** 2025–2026

---

## 📝 License

This project is released under the MIT License — feel free to use, modify, and distribute with attribution.

---

## 🤝 Acknowledgments

Built upon foundational research in joint radar-communication systems by leading researchers in the ISAC field. Special thanks to the IEEE Signal Processing community for open-access publications that made this work possible.

---

⭐ **If you find this project useful, please star the repository!**
