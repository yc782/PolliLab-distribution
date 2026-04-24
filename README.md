# PolliLab
This script is intended to be used with tablet image data obtained from ImageJ processing, where the csv dataframe consists of particle number, their respective area (micron), XY coordinates (centroids) and mean, min, max (not useful information. those are stats of pixel color). The first row in the dataframe contains the are and XY coordinate of the tablet (circle ROI), and the last four rows are the summary output from ImageJ. The last 4 rows will be ignored.

This program divide the tablet into N sectors of equal area (16 by default), counts the number of particles in each sector, summarizes area and diameter of each sector, and provides particle size distribution of the tablet (global) in 10 bins by default unless specified by max and mix values or number of size bins. 
The script produces .xlsx file that contains: 
1) Summary sheet with sector stats, global area bin counts, global diameter bin counts, D10/D50/D90 using diamater, log-normal fit using diamater
2) Sector sheets with raw particle data, sector summary, global area distribution, local area distribution, global diameter distribution, local diameter distribution and histograms for each distribution.

Note: 
- Global bins use the same bin edges across all sectors.
- Local bins are re-created separately inside each sector.
- With --log_bins, both global and local bins use logarithmic spacing
- Log-normal fitting uses only positive area values

```
$ python distribution.py -h
usage: distribution.py [-h] [-o O] [--n_bin N_BIN] [--max_area MAX_AREA] [--min_area MIN_AREA] [--min_diameter MIN_DIA] [--max_diameter MAX_DIA] [--n_bin_size N_BIN_SIZE] [--log_bins] input

Compute statistics about particle distribution

positional arguments:
  input                 Input file, csv format

options:
  -h, --help            show this help message and exit
  -o O                  Optional name of output file, csv format
  --n_bin N_BIN         Number of Sector-shaped bins; default=16
  --max_area MAX_AREA   Optional max Area for the size distribution (auto-calculated if not supplied)
  --n_bin_size N_BIN_size    Number of size bins for area/diameter; default=10
  --max_diameter MAX_DIA    Optional maximum diameter for global diameter bins
  --log_bins            Use logarithmic bins for area and diameter distributions
