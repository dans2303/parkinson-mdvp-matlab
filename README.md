# 🧠 Parkinson Detection using MDVP Features (MATLAB)

## 📌 Project Overview

This project analyzes voice recordings using medically relevant MDVP (Measures of Dysphonia in Voice Profile) features to support the early detection of Parkinson’s Disease (PD). The work is structured in two parts:

1. **MATLAB-based clinical signal analysis and feature engineering**
2. **(Coming Soon)** Machine Learning-based classification using MATLAB and Python

The goal is to extract and interpret acoustic biomarkers that are sensitive to neuromuscular degradation, such as pitch and amplitude instability.

---

## 🧮 Part I – MATLAB-Based Voice Analysis & Feature Engineering

### 🎯 Objective

To detect PD symptoms by analyzing variation in voice features such as jitter, shimmer, and harmonic-to-noise ratio. MATLAB is used to perform statistical calculations and generate diagnostic plots.

### 📊 Dataset Summary

- 198 `.wav` voice recordings  
- Extracted features include:
  - **Pitch**: `AVG_F0`, `STD_F0`, `PFR`
  - **Jitter**: `JITA`, `JITT`, `RAP`, `PPQ`
  - **Shimmer**: `SHIM`, `SHDB`, `APQ`, `SAPQ`
  - **Noise Measures**: `NHR`, `VTI`, `VAM`, `NSH`, `DSH`
- `Status`: Binary label (0 = Healthy, 1 = Parkinson's)

### 📘 Theory Behind MDVP Features

To understand the medical and signal processing basis behind the features used, see the [MDVP Feature Theory](mdvp_theory.md).

---

## 🤖 Part II – Machine Learning Classification (Coming Soon)

This section will use engineered features to train classifiers such as:
- XGBoost
- Support Vector Machine (SVM)
- Linear Discriminant Analysis (LDA)
- K-Nearest Neighbors (KNN)

Evaluation will include accuracy, precision-recall, and confusion matrices.

---

## 📁 Project Structure

<pre><code> ``` parkinson-mdvp-matlab/ ├── code/ # MATLAB scripts for feature analysis ├── data/ # Voice dataset (CSV or MAT) ├── results/ # Plots, correlation matrices, PCA, etc. ├── mdvp_theory.md # In-depth explanation of MDVP features ├── README.md # This documentation ``` </code></pre>


---

## ▶️ How to Run (MATLAB)

1. Open `M1.m` in MATLAB
2. Make sure `data/` contains the voice dataset
3. Run the script to generate plots and statistical summaries in `/results`

---

## 📄 License

This project is licensed under the **MIT License**.  
Feel free to use, adapt, and reference with proper attribution.

---

## 📫 Contact

Mirna Danisa Tandjung  
Email: [mirnadanisat@gmail.com]  
LinkedIn: [(https://www.linkedin.com/in/mirnadanisatandjung/)]
