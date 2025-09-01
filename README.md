# Credit Counterparty Risk in the LGM Model

This project explores **counterparty credit risk (CCR)** in interest rate derivatives by comparing two one-factor short-rate models: **Hull–White (HW)** and **Linear Gaussian Markov (LGM)**.  

---

## Overview

We study three **portfolios of interest rate swaps** and investigate how CCR metrics evolve under different modeling and collateralization assumptions.  

The project also provides an **analytical extension of the Hull–White model** and a **derivation of the main properties of the LGM model**, showing how it generalizes HW while keeping a tractable structure for calibration and pricing.

Both models are calibrated to the same swaption data and used in **Monte Carlo simulations** to produce exposure profiles.

---

## Methodology

1. **Model Foundations**
   - Hull–White calibrated with mean reversion `a` and constant volatility `σ`.
   - LGM derived as an extension of HW, with deterministic (piecewise-constant) volatility.
   - Analytical derivation of bond prices and swaption valuation via Jamshidian’s decomposition.

2. **Calibration**
   - HW: joint optimization of `(a, σ)`.
   - LGM: bootstrap calibration of piecewise-constant volatilities.
   
3. **Exposure Simulation**
   - Monte Carlo horizon: **10 years**.
   - Time grids: **quarterly**, **weekly**, **daily**.
   - Portfolios simulated both **with and without netting/collateral**.

4. **CCR Metrics**
   - **Expected Exposure (EE)**
   - **Expected Positive Exposure (EPE)**
   - **Potential Future Exposure (PFE)**
   - **Peak PFE**

---

## Results

- **Model Comparison**  
  Deviations between HW and LGM are **below 5%** across all horizons, except for the very short term.  
  The added flexibility of LGM does **not materially change CCR profiles** in this setup.

- **Collateral Effects**
  - **Semiannual collateral** reduces exposures but does not eliminate risk.  
  - **Weekly collateral** drives exposures almost to zero.  
  - **Daily grid with weekly collateral** increases runtime but remains consistent.

- **Portfolios**  
  Netting significantly lowers exposures. Without netting, each trade is treated as standalone, inflating risk.

---

## Key Insights

- **Hull–White vs LGM**  
  For standard swap portfolios and smooth volatility surfaces, HW is nearly as effective as LGM despite being simpler.  

- **Analytical Value of LGM**  
  The LGM framework extends HW analytically, preserving tractability while allowing for **time-dependent volatility**.  
  This makes it a stronger candidate in **irregular markets** or for **optional/path-dependent instruments**.  

- **Complexity vs Computation**  
  LGM’s extra complexity is **conceptual, not computational**.  

---

## Implementation Notes

- Implementation in **MATLAB**.  
- Monte Carlo with up to **250,000 scenarios**.  
- Optimizations:
  - Affine parameter reuse  
  - Selective numerical integration  
  - Adjustments for finer grids  

---

  
