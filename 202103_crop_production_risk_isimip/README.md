# Assessing crop production risks under 21st century warming


This repository contains a folder for the code used in a research article. The code replicates the figures and tables, and provides calibration code implemented in [CLIMADA Python](https://github.com/CLIMADA-project/climada_python) using input data from [ISIMIP](https://www.isimip.org/) and [FAOSTAT](http://www.fao.org/faostat/).


*Paper in preparation.*
Corresponding author: Samuel Eberenz


## Most important scripts:

### isimip3b_crop_config.py
(python script)

configure parameters, directories, variables, etc.

### main_isimip3b.py
(python script)

Wrapper script, run to initiate CLIMADA hazard and exposure sets, compute impacts, bin impacts by global mean temperature (GMT) and calculate and export statistics.
Results, among intermediate outputs, in statistics per country and crop type; and main results table as CSV

*Requires:*
* isimip3b_crop_config as co
* isimip3b_climada_wrappers as cw (--> requires CLIMADA v1.5.1+)
* isimip3b_gmt_binning as gmt_binning
* isimip3b_impact_statistics as impact_statistics

### make_crop_production_risk_world_map_plots.ipynb
(jupyter notebook)

Make world maps of crop production statistics, both gridded and per country.

*Requires:*
* isimip3b_crop_config
* plot_worldmaps_utils.py
* cartopy
* mplotutils (https://github.com/mathause/mplotutils)

### make_country_stat_plots.py
(python script)

Make plots per country for production deviation per return periods (RP), Probability ratio (PR), probability (PP), Coefficient of Variation (CV).
*Requires:*
* isimip3b_crop_config
* plot_stats_per_country_utils

## Notes:
Please note: most contents of the folders 'data' and 'results' is currently not synchronized. Change '.gitignore' to track the content of these folders.
