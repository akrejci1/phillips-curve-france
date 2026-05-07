# New Keynesian Phillips Curve – France

Econometrics assignment estimating the New Keynesian Phillips Curve (NKPC) for France using quarterly OECD data from 1980 Q1 to 2019 Q4. The analysis replicates and extends the framework of Galí & Gertler (1999) using GMM, NLS, and ML estimation methods.

---

## Overview

The NKPC links current inflation to expected future inflation and a measure of real marginal costs. This project estimates both reduced-form and structural versions of the curve for France, compares ULC and output gap as proxies for real marginal costs, tests for hybrid (backward-looking) dynamics, and assesses long-run vertical Phillips curve restrictions.

The sample ends at 2019 Q4 to avoid distortions caused by the COVID-19 pandemic.

---

## Data

All data are sourced from [OECD Data Explorer](https://data-explorer.oecd.org/).

| Variable | Description | Unit |
|---|---|---|
| `gdp_real` | Real GDP (chain-linked) | Mil. EUR, 2020=100 |
| `gdp_nom` | Nominal GDP | Mil. EUR, current prices |
| `ulc` | Unit Labour Cost index | Index, 2015=100 |
| `ir_10y` | Long-term interest rate (10Y) | % p.a. |
| `ir_3m` | Short-term interest rate (3M) | % p.a. |

**Required files** (place in the working directory):
- `gdp.xlsx` – sheet `Table`, row 1: dates, row 2: real GDP, row 3: nominal GDP
- `interest_rates.xlsx` – sheet `Table`, row 2: 10Y rate, row 3: 3M rate
- `unit_labour_cost.xlsx` – column 2: ULC index values

---

## Model Variables

| Variable | Construction |
|---|---|
| `pi_def` | Quarterly inflation: `diff(log(GDP deflator))` |
| `gap` | Output gap: HP-filter cycle on `log(gdp_real)`, λ=1600 |
| `ulc_c` | Cyclical ULC: HP-filter cycle on `log(ulc)`, λ=1600 — proxy for real marginal costs |
| `spread` | Interest rate spread: `(ir_10y − ir_3m) / 100` |
| `pi_w` | Wage inflation: `diff(log(ulc))` |

---

## Models

All GMM models use lags 1–4 of all five model variables as instruments (20 instruments total) and a HAC covariance matrix (Quadratic Spectral kernel).

### Reduced-Form NKPC

$$\pi_t = \lambda x_t + \beta \mathbb{E}_t[\pi_{t+1}] + u_t$$

| Model | Regressor | λ | β | J-test p |
|---|---|---|---|---|
| Model 1 | ULC cycle | 0.044** | 1.000 | 0.160 |
| Model 2 | Output gap | 0.011 | 0.999 | 0.038 |

Model 1 (ULC) passes the J-test and yields a positive, significant λ. Model 2 (gap) fails the J-test, indicating misspecification.

### Structural NKPC (Calvo pricing)

$$\pi_t = \lambda mc_t + \beta \mathbb{E}_t[\pi_{t+1}]$$
$$\lambda = \frac{(1-\theta)(1-\beta\theta)}{\theta}$$

Two normalisations of the moment conditions are estimated (eqs. 18 and 19 in GG1999):

| Model | Normalisation | θ | β | λ | Fixation (Q) | J-test p |
|---|---|---|---|---|---|---|
| Model 3 | (18) | 0.779 | 1.005 | 0.061 | 4.53 | 0.326 |
| Model 4 | (19) | 0.811 | 1.000 | 0.044 | 5.28 | 0.160 |
| GG1999 | (18) | 0.829 | 0.926 | 0.047 | 5.85 | — |
| GG1999 | (19) | 0.884 | 0.941 | 0.021 | 8.62 | — |

### Hybrid NKPC

$$\pi_t = \lambda mc_t + \gamma_f \mathbb{E}_t[\pi_{t+1}] + \gamma_b \pi_{t-1} + u_t$$

| Model | Form | λ | γ_f | γ_b | J-test p |
|---|---|---|---|---|---|
| Model 5 | Reduced | −0.005 | 0.586 | 0.402 | 0.001 |
| Model 6 | Structural | — | — | — | — |

The ratio γ_f/γ_b ≈ 1.4 suggests forward-looking behaviour dominates. However, both hybrid models fail the J-test.

### Restricted Models (β = 1, long-run verticality)

| Model | θ | λ | J-test p |
|---|---|---|---|
| Model 4 unrestricted (β ≈ 0.997) | 0.810 | 0.044 | 0.160 |
| Model 4b restricted (β = 1) | 0.806 | 0.047 | 0.183 |

Imposing β = 1 barely changes structural estimates, consistent with the unrestricted β already being ≈ 1.

---

## Alternative Estimators

| Method | θ | β | λ | Fixation (Q) | Log-lik |
|---|---|---|---|---|---|
| GMM (Model 4) | 0.811 | 1.000 | 0.044 | 5.28 | — |
| NLS | 0.788 | 0.931 | 0.071 | 4.72 | 656.3 |
| ML (normal) | 0.788 | 0.931 | 0.071 | 4.73 | 656.3 |
| ML (t-dist.) | 0.808 | 0.930 | 0.059 | 5.22 | 659.1 |

NLS and ML produce identical results (as expected). The systematically lower β ≈ 0.931 relative to GMM reflects their inability to handle the endogeneity of expected inflation.

### LR Test for Long-Run Verticality (H₀: β = 1)

| Method | LR stat | p-value | Decision (5%) |
|---|---|---|---|
| NLS | 4.530 | 0.033 | Reject H₀ |
| ML (normal) | 4.530 | 0.033 | Reject H₀ |
| ML (t-dist.) | 3.909 | 0.048 | Reject H₀ |

All three methods narrowly reject H₀ at 5% — but this is an artefact of endogeneity bias pulling β below 1. GMM, which handles endogeneity, leaves the restriction essentially binding.

---

## Diagnostics

Residual autocorrelation is checked via Ljung-Box tests and ACF plots.

| Model | LB(8) p | LB(12) p |
|---|---|---|
| Model 3 (full sample) | 0.002 | 0.006 |
| Model 4 (full sample) | 0.001 | 0.005 |
| Model 5 (hybrid) | <0.001 | <0.001 |
| Model 3 (1980–2000) | 0.078 | 0.202 |

The narrow pre-euro sub-sample (1980–2000) eliminates residual autocorrelation in Model 3. The best overall compromise for the full sample is **Model 4**, where the J-test passes and the HAC covariance matrix robustly handles remaining autocorrelation.

---

## Key Findings

- **ULC as the main driver of inflation:** The cyclical component of unit labour costs is the only statistically significant proxy for real marginal costs in France (λ = 0.044, p = 0.012). The output gap fails both on significance and instrument validity.
- **Price stickiness:** Structural estimates imply firms keep prices fixed for approximately 5 quarters on average, in line with the GG1999 benchmark.
- **Endogeneity matters:** The disparity between GMM (β ≈ 1.000) and NLS/ML (β ≈ 0.931) is a textbook illustration of the bias from ignoring endogeneity of expected inflation.
- **Heavy tails:** An ML model with Student-t errors significantly outperforms the normal specification (LR = 5.53, p = 0.019), with estimated degrees of freedom ν ≈ 5.36.
- **Residual autocorrelation:** Persistent autocorrelation in structural NKPC residuals likely reflects structural features of the French labour market and the transition to the eurozone.


---

## Requirements

```r
install.packages(c("openxlsx", "ggplot2", "ggfortify", "mFilter",
                   "forecast", "gmm", "minpack.lm", "maxLik"))
```

---

## Reference

Galí, J., & Gertler, M. (1999). Inflation dynamics: A structural econometric analysis. *Journal of Monetary Economics*, 44(2), 195–222.
