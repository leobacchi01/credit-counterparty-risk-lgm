# credit-counterparty-risk-lgm
Monte Carlo simulation of Credit Counterparty Risk under LGM and Extended Vasicek model.




This work analyzes counterparty credit risk (CCR) for three portfolios of interest
rate swaps using two one-factor short-rate models: Hull–White (HW) and Linear
Gaussian Markov (LGM). Both models are calibrated to the same set of at-the-money
swaptions.
In our implementation, the HW model jointly calibrates the mean reversion speed
a and a single volatility level, while the LGM uses a piecewise-constant volatility
structure, keeping a fixed from the outset.
Monte Carlo simulations were run over a 10-year horizon with different time grids:
a quarterly grid to model collateral re-margining on a semi-annual basis, and weekly
and daily grids to simulate collateral paid weekly. These simulations produce
the main exposure profiles — expected exposure (EE), expected positive exposure
(EPE), potential future exposure (PFE) and peak PFE — both with and without
netting and collateral agreements.
Across all configurations, the two models show deviations below 5%, with a more
noticeable difference only at the first time horizon.
The results indicate that, for standard swap portfolios and a relatively smooth
volatility surface as observed in the available data, the greater structural flexibility
of the LGM does not lead to significant differences in counterparty risk metrics.
However, the more detailed parameterization does not result in higher calibration
or simulation times: the additional burden is conceptual, not computational.
The LGM remains a more general extension of the HW model and may offer advantages
in markets characterized by strong curve irregularities or in the valuation of
optional or path-dependent instruments.
