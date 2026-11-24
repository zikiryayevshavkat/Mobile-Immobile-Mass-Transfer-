# Suspension Filtration Reproduction

This archive contains the model, data, and figure scripts for Figures 1–10 and Tables 1–7 of the revised manuscript. Figure 1 is a schematic. Plot labels use words, with the corresponding mathematical notation and units given in the manuscript captions.

Run these commands from the extracted directory with Python 3.11 or later. The numerical results were generated with Python 3.12.14 and NumPy 2.3.5.

```text
python -m pip install -r requirements.txt
python solver.py
python plot_results.py
python sensitivity.py
python plot_sensitivity.py
python digitize_column_observations.py
python column_validation.py
python plot_column_validation.py
python column_straining.py
python plot_straining.py
```

The solver defines the published parameters in `Parameters`. The examples use a 0.30 mobile porosity and a 0.20 immobile porosity, a fixed Darcy velocity, and a 0.007 deposit threshold. These are illustrative parameter choices. The rate sweep holds these parameters fixed while varying only the exchange coefficient.

The solver applies explicit exchange and deposition updates, credits or debits the mobile concentration by the actual storage changes, and integrates the transport fluxes implicitly. The inlet has a fixed face concentration and the outlet carries only the advective flux.

The script reruns all four exchange cases, separate spatial and temporal refinement studies, an outlet-distance check, and the Ogata–Banks benchmark with exchange and deposition disabled. The benchmark verifies the transport operator. It uses analytical cell averages from three-point quadrature.

- `results/profiles.csv` contains the cell coordinates and the three concentrations for every exchange case and saved time.
- `results/profiles.npz` contains the same profiles and the numerical and analytical benchmark curves for plotting. Its `fields` axes are exchange case, time, concentration field, and spatial cell.
- `results/verification.json` contains the parameter comparison, mass budgets, refinement changes, and benchmark errors. Relative values are stored as fractions.
- `results/sensitivity_budgets.csv` gives the integrated budgets and deposit centroids for the 16 production cases at the four saved times.
- `results/sensitivity.json` includes those budgets, point concentrations, and the four additional runs on the finer grid. Case keys are ordered as immobile pore fraction, exchange rate, and deposition threshold.
- `results/sensitivity_profiles.npz` gives the spatial fields for the sensitivity production runs in case, time, concentration field, and spatial cell order.
- `results/column_observations.csv` gives the plotted experimental-marker centers from Wang et al. (2014), DOI `10.1016/j.crme.2014.07.001`, for Figs. 2, 3a, 3b, 5a, and 6a. `digitize_column_observations.py` records the page, rendering scale, crop and coordinate transform. The original figures and their fitted curves are not redistributed.
- `results/column_validation.json` gives column properties, parameter estimates, search bounds, model curves and error metrics at both resolutions. `results/column_predictions.csv` matches every observed value to each model prediction. The 0.25 g/L outlet points through 60 pore volumes are used in fitting. Later outlet points at that feed, all the deposit points, and all retained 0.50 g/L outlet points are held out.
- `column_straining.py` adds a depth factor to both capture branches, with fixed median grain diameter of 0.000410 m and attenuation exponent of 0.432. The solver defaults to uniform capture. The original fit split is retained and the early rate fit is interpolated across capture fractions. The fully depth-dependent capture case fixes the inlet straining fraction to one as a structural diagnostic with prescribed exchange.
- `results/straining_validation.json` gives the fitted parameters, refinement checks and metrics for the revised models. `results/straining_predictions.csv` pairs every digitized observation with each revised model's prediction. The deposit residuals motivated the extension, and every deposit point stayed outside numerical calibration. The comparison is retrospective.
- The ten figures are available in PNG and vector PDF formats.
