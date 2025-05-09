# STABLE: SubTropical Atmospheric ridge and BLocking Events
(STABLE has been renamed from a previous version, from BLOCS)

An open-source and user-friendly Python algorithm for detecting and tracking atmospheric blockings and subtropical ridge obstruction events.
This is a working copy in python 3.9 of the methodology presented in [Sousa et al. (2021)](https://doi.org/10.1175/JCLI-D-20-0658.1). Some alterations are available as a namelist_input.txt file.
Questions, suggestions, and corrections: Miguel M. Lima (malima@ciencias.ulisboa.pt).

[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# What do I need to get and run STABLE?

## To run STABLE, you need

   * [![python](https://img.shields.io/badge/Python-3-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)
   * [![Git: ](https://img.shields.io/badge/Git--blue)](https://git-scm.com/)

and

 * [![Anaconda 3](https://img.shields.io/badge/Anaconda-3-green.svg)](https://www.anaconda.com/) (or similar to manage python packages)

or

  *  [![python](https://img.shields.io/badge/Python-3-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org) and the required modules on a cluster

## The packages required to run STABLE are:
  
```
- numpy
- xarray
- scipy
- os
- pickle
- pandas
- tqdm
```

# Installation

1 - Clone STABLE repository.

 ```
git clone https://github.com/mikaslima/STABLE.git
  ```

2 - Verify you have installed all packages requiered for STABLE. If you use an Anaconda environment, please be sure you have activated the environment.

# Running the algorithm

## Input data

The code takes any resolution, regularly spaced data (tested with 2.5º, 1º, and 0.25º). Ideally, the name should be "Z500_(start year)_(end year)_(hemisphere)_(name of institution)" for example: "Z500_1940_2023_NH_ECMWF.nc". The variable names inside the file should be "z" for the geopotential at 500 hPa in units of [m], "time" as datetime format [YYYY-MM-DD], "lat", and "lon" for the coordinates in degrees.

In this release, you can find a 2-year (2019-2020) dataset from NCAR, at 2.5º resolution, for the Northern Hemisphere, ready to run. This example data returns a catalogue that can be seen in the ```Data/Output_data```. Additionally, example figures are shown in ```Figures```. Larger 80-year datasets are available at [![Zenodo: https://doi.org/10.5281/zenodo.13891996](https://img.shields.io/badge/Zenodo-10.5281/zenodo.13891996-blue)](https://doi.org/10.5281/zenodo.13891996).

## Code

The ```Algorithm``` itself is divided in 3 major steps/scripts: Daily structure identification, structure tracking in time, and production of the catalogues with general statistic (by daily observation, event) and the respective observation masks.
In each of these, the preamble contains a set of variables to be changed by the user in the "Data/Input_data/namelist_input.txt" file (see description below).
These codes may have significant departures from the MATLAB application of [Sousa et al. (2021)](https://doi.org/10.1175/JCLI-D-20-0658.1), the ones I considered are marked as such in the code.
To run successfully the code please edit the namelist according to your needs and run the scripts in succession.

# Namelist description

Name                |  Type         |   Values    |   Replicate Sousa et al. (2021)      |     Description
--------------------|---------------|-------------|--------------------------------------|----------------------
year_file_i         |  int          |   any       |   any                                |    first year of the data file
year_file_f         |  int          |   any       |   any                                |    last year of the data file
date_init           |  str          |   any       |   '1950-01-01'                       |    start date of the analysis
date_end            |  str          |   any       |   '2020-12-31'                       |    end date of the analysis
res                 |  float        |   any       |   2.5                                |    resolution of the data (e.g., 2.5, 1, 0.25)
region              |  string       |   NH, SH    |   NH, SH                             |    (NH) Northern Hemisphere<br>(SH) Southern Hemisphere
data_type           |  string       |   any       |   NCAR                               |    name of institution of the data, in the name of the file
min_struct_area     |  int          |   any       |   500000                             |    threshold total area for the structures (usually 500000)
use_max_area        |  int          |   0,1       |   0                                  |    (0) Don't use a maximum area threshold (1) Use a maximum area threshold
max_struct_area     |  int          |   any       |   -                                  |    maximum threshold total area for the structures
n_days_before       |  int          |   any       |   15                                 |    number of days before the analysis day to compute the LATmin (usually 15)
assymetrical_LATmin |  int          |   1,2       |   1                                  |    (1) Uses a horizontal LATmin as in Sousa et al. (2021)<br>(2) Uses a variable LATmin with longitude
omega_hybrid_method |  int          |   1,2       |   1                                  |    (1) Uses the same condition to identify hybrid blocks as in Sousa et al. (2021)<br>(2) Uses a new condition, which considers the Rex area within the structure.
GHGN_condition      |  int          |   1,2,3     |   1                                  |    (1) Considers polar area the same way as in Sousa et al. (2021)<br>(2) Computes the gradient northward of the (90-delta) boundary and extends the structure conditions over this boundary, but cut the area poleward of 85º for continuity problems<br>(3) Completely ignores the area poleward of (90-delta)
lat_polar_circle    |  float or int |   any       |   75                                 |    Latitude of the polar circle to consider a Rex block to be of a polar nature (ideally between ~60 up to 90, usually is set to 75)
delta               |  float or int |   any       |   15                                 |    Delta in degrees to compute the gradients (usually set to 15 degrees)
tracking_method     |  int          |   1,2       |   1                                  |    (1) Use a "OR" condition, either the overlapping area is exceeding the "area_threshold" in the day of the analysis or the previous<br>(2) Use a "AND" condition, the overlapping area must exceed the "area_threshold" in both of the day of the analysis and the previous.
overlap_threshold   |  float        |   ]0,1[     |   0.5                                |    Fraction of overlap needed to consider the evolution of the structure (usually set at 0.5)
persistence         |  int          |   any       |   4                                  |    Minimum number of days for an atmospheric blocking to be considered an event
full_ridges         |  int          |   0,1       |   0                                  |    (0) Consider events that only have subtropical ridges as observations (1) Disregard these types of events
full_polar          |  int          |   0,1       |   0                                  |    (0) Consider events that only have polar blocks (1) Disregard these types of events
catalogue_output    |  int          |   0,1       |   0                                  |    (0) Catalogue full events and their daily components (1) Catalogue single occurrences only
get_masks           |  int          |   0,1       |   0,1                                |    (0) Disregard masks (1) Retrieve the masks of the events and the 2D intensity index
save_latmin         |  int          |   0,1       |   0,1                                |    (0) Do not save LATmin data (1) Save LATmin data as netcdf
get_type            |  int          |   0,1       |   0,1                                |    (0) Do not produce netcdf with blocking types (1) Save netcdf with local, structure, and event type: the hundreth represents the event type (100-Ridge; 200-Omega; 300-Hybrid; 400-Rex; 500-Polar), the tenth the daily structure type (10-Ridge; 20-Omega; 30-Hybrid; 40-Rex; 50-Polar), and the unit is the local type (1-Ridge; 2-Omega; 3-Rex)
-----------------------------------------------------------------------------------------------------------------

# Post processing examples

Also available are a series of figures to study the obtainable catalogues, in ```Post_processing```. For example, the climatology obtained through the given example dataset:

![plot](./Figures/02-Rec_Fig7(climatology).jpg)
