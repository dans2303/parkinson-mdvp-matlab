# 🧠 MDVP Feature Theory for Parkinson’s Detection

This document provides a theoretical background on the key acoustic features extracted from voice recordings using MDVP
---

## 🎯 Objective

To extract and interpret medically-relevant voice features that are sensitive to neuromuscular changes, enabling classification of Parkinson’s patients from healthy controls.

The following features are the most significant for Parkinson’s Disease detection:

---

## 🔹 1. Jitter (JITT, JITA, RAP, PPQ)

Jitter measures the **frequency instability** of vocal fold vibrations — it reflects how much the pitch (Fo) varies between cycles.

- **Why it matters**: Parkinson’s Disease often causes vocal tremor, which leads to irregular pitch.
- **Used Features**: `JITT`, `JITA`, `RAP`, `PPQ`, `SPPQ`

**Formula (Jitter %):**
\[
Jitter(\%) = \frac{1}{N-1} \sum_{i=1}^{N-1} \left| \frac{T_i - T_{i+1}}{T_i} \right| \times 100
\]
Where \( T_i \) is the duration of the \( i^{th} \) pitch period.

---

## 🔹 2. Shimmer (SHIM, SHDB, APQ)

Shimmer captures **amplitude instability** — fluctuations in the loudness of the voice across vocal cycles.

- **Why it matters**: Amplitude irregularity reflects vocal fatigue or instability, common in PD.
- **Used Features**: `SHIM`, `SHDB`, `APQ`, `SAPQ`

**Formula (Shimmer %):**
\[
Shimmer(\%) = \frac{1}{N-1} \sum_{i=1}^{N-1} \left| \frac{A_i - A_{i+1}}{A_i} \right| \times 100
\]
Where \( A_i \) is the peak amplitude of the \( i^{th} \) cycle.

---

## 🔹 3. Harmonics-to-Noise Ratio (NHR)

NHR quantifies the ratio of periodic (harmonic) content to aperiodic (noisy) content in the voice.

- **Why it matters**: Lower NHR indicates more breathiness and vocal noise — both are signs of impaired vocal cord closure in PD.

**Formula:**
\[
\text{NHR} = \frac{\text{Harmonic Energy}}{\text{Noise Energy}}
\]

---

## 🔹 4. Voice Tremor Indicator (VTI)

VTI reflects **low-frequency modulation** in amplitude — a proxy for vocal tremor.

- **Why it matters**: Parkinson’s tremors can manifest in sustained vowel sounds as periodic volume fluctuations.

---

## 🔹 5. Standard Deviation of Fundamental Frequency (STD_F0)

STD_F0 shows how much the fundamental frequency varies over time.

- **Why it matters**: Healthy individuals have smoother Fo variation, while PD patients often display irregular pitch control.

---

## 🔹 6. Pitch Frequency Range (PFR)

PFR is the difference between the maximum and minimum pitch during sustained phonation.

- **Why it matters**: Monotone speech or abnormal pitch variation is a known symptom of PD.

---

## 🗂️ Summary Table of Feature Categories

| Category        | Features Included                   | Clinical Relevance                       |
|-----------------|-------------------------------------|------------------------------------------|
| Jitter          | JITT, JITA, RAP, PPQ, SPPQ          | Frequency instability, vocal tremor      |
| Shimmer         | SHIM, SHDB, APQ, SAPQ               | Amplitude instability                    |
| Noise Ratio     | NHR, NSH, DSH                       | Voice breathiness and glottal noise      |
| Frequency Stats | AVG_F0, STD_F0, PFR, VF0            | Pitch control and variability            |
| Modulation      | VAM, VTI, DVB                       | Tremor-like amplitude changes            |

---

📌 These features were extracted from `.wav` recordings using MATLAB scripts. They are later used for classification models in the machine learning phase of the project.
