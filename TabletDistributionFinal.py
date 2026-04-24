# This script produce:
#   Summary sheet: 1.sector stats; 2. global Area-bin counts; 3. global Diameter-bin counts
#                 4. D10/D50/D90 using diameter; 5. log-normal fit using Diameter
#   Each sector sheet: 1. raw particle data; 2. sector summary; 3. global Area-bin distribution
#                    4. local Area-bin distribution; 5. global Diameter-bin distribution
#                    6. local Diameter-bin distribution; 7. histograms for each distribution

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Compute particle distribution statistics by sector')
parser.add_argument('input', help='Input file, csv format')
parser.add_argument('-o', help='Optional output filename, xlsx format')
parser.add_argument('--n_bin', type=int, default=16, help='Number of angular sector bins; default=%(default)s')
parser.add_argument('--n_bin_size', type=int, default=10, help='Number of size bins for area/diameter; default=%(default)s')
parser.add_argument('--min_area', type=float, default=None, help='Optional minimum area for global area bins')
parser.add_argument('--max_area', type=float, default=None, help='Optional maximum area for global area bins')
parser.add_argument('--min_diameter', type=float, default=None, help='Optional minimum diameter for global diameter bins')
parser.add_argument('--max_diameter', type=float, default=None, help='Optional maximum diameter for global diameter bins')
parser.add_argument('--log_bins', action='store_true', help='Use logarithmic bins for area and diameter distributions')
args = parser.parse_args()


def compute_percentiles(series):
    s = pd.to_numeric(series, errors='coerce').dropna()
    if len(s) == 0:
        return np.nan, np.nan, np.nan
    return np.percentile(s, 10), np.percentile(s, 50), np.percentile(s, 90)


def lognormal_fit(series):
    s = pd.to_numeric(series, errors='coerce').dropna()
    s = s[s > 0]
    if len(s) < 2:
        return np.nan, np.nan
    log_s = np.log(s)
    mu = log_s.mean()
    sigma = log_s.std(ddof=1)
    return mu, sigma


def lognormal_pdf(x, mu, sigma):
    x = np.asarray(x)
    pdf = np.full_like(x, np.nan, dtype=float)
    if np.isfinite(mu) and np.isfinite(sigma) and sigma > 0:
        valid = x > 0
        xv = x[valid]
        pdf[valid] = (1.0 / (xv * sigma * np.sqrt(2 * np.pi))) * np.exp(-((np.log(xv) - mu) ** 2) / (2 * sigma ** 2))
    return pdf


def make_bins(data, n_bins, use_log=False, min_val=None, max_val=None):
    s = pd.to_numeric(data, errors='coerce').dropna()
    if use_log:
        s = s[s > 0]

    if len(s) == 0:
        return None

    lo = s.min() if min_val is None else min_val
    hi = s.max() if max_val is None else max_val

    if use_log:
        if lo <= 0:
            positive = s[s > 0]
            if len(positive) == 0:
                return None
            lo = positive.min()
        if hi <= 0 or lo >= hi:
            return None
        return np.logspace(np.log10(lo), np.log10(hi), n_bins + 1)
    else:
        if not np.isfinite(lo) or not np.isfinite(hi):
            return None
        if lo == hi:
            hi = lo + 1e-9
        return np.linspace(lo, hi, n_bins + 1)


def make_distribution_table(series, bins, label):
    if bins is None or len(series.dropna()) == 0:
        return pd.DataFrame(columns=[label, 'Count'])

    cut_series = pd.cut(series, bins=bins, include_lowest=True)
    counts = cut_series.value_counts().sort_index()
    out = counts.reset_index()
    out.columns = [label, 'Count']
    return out


# -------------------------
# Validate input
# -------------------------
if not args.input.lower().endswith('.csv'):
    raise SystemExit('ERROR: Input must be a CSV file.')

# -------------------------
# Load data
# -------------------------
df = pd.read_csv(args.input, skipfooter=4, engine='python')

metadata = df.iloc[0]
df = df.iloc[1:].copy()

for col in ['X', 'Y', 'Area']:
    df[col] = pd.to_numeric(df[col], errors='coerce')

df = df.dropna(subset=['X', 'Y', 'Area']).copy()

x0 = float(metadata['X'])
y0 = float(metadata['Y'])

# Shift coordinates
df['X'] = df['X'] - x0
df['Y'] = df['Y'] - y0

# Compute diameter from area
df['Diameter'] = np.sqrt(4.0 * df['Area'] / np.pi)

# Polar coordinates
df['r'] = np.sqrt(df['X'] ** 2 + df['Y'] ** 2)
df['theta'] = np.arctan2(df['Y'], df['X'])

# Sector assignment
n_bin = args.n_bin
bin_col = f'bin{n_bin}'
df[bin_col] = np.floor((df['theta'] + np.pi) / (2 * np.pi) * n_bin) % n_bin
df[bin_col] = df[bin_col].astype(int)

# -------------------------
# Global bins
# -------------------------
area_source = df['Area'].copy()
diam_source = df['Diameter'].copy()

if args.log_bins:
    area_source = area_source[area_source > 0]
    diam_source = diam_source[diam_source > 0]

global_area_bins = make_bins(
    area_source,
    args.n_bin_size,
    use_log=args.log_bins,
    min_val=args.min_area,
    max_val=args.max_area
)

global_diameter_bins = make_bins(
    diam_source,
    args.n_bin_size,
    use_log=args.log_bins,
    min_val=args.min_diameter,
    max_val=args.max_diameter
)

if global_area_bins is None:
    raise SystemExit('ERROR: Could not create global area bins.')
if global_diameter_bins is None:
    raise SystemExit('ERROR: Could not create global diameter bins.')

df_area = df.copy()
df_diam = df.copy()

if args.log_bins:
    df_area = df_area[df_area['Area'] > 0].copy()
    df_diam = df_diam[df_diam['Diameter'] > 0].copy()

df_area['global_area_bin'] = pd.cut(df_area['Area'], bins=global_area_bins, include_lowest=True)
df_diam['global_diameter_bin'] = pd.cut(df_diam['Diameter'], bins=global_diameter_bins, include_lowest=True)

# -------------------------
# Summary stats by sector
# -------------------------
stats = df.groupby(bin_col).agg(
    particle_count=('Area', 'count'),
    sum_area=('Area', 'sum'),
    mean_area=('Area', 'mean'),
    std_area=('Area', 'std'),
    sem_area=('Area', 'sem'),
    mean_diameter=('Diameter', 'mean'),
    std_diameter=('Diameter', 'std'),
    sem_diameter=('Diameter', 'sem')
)

sector_percentiles = df.groupby(bin_col)['Diameter'].apply(
    lambda s: pd.Series(compute_percentiles(s), index=['D10_diameter', 'D50_diameter', 'D90_diameter'])
).unstack()

stats = stats.join(sector_percentiles)

global_area_distr = (
    df_area.groupby([bin_col, 'global_area_bin'], observed=False)
    .size()
    .unstack(fill_value=0)
)
global_area_distr.columns = [f'Area_{c}' for c in global_area_distr.columns]

global_diameter_distr = (
    df_diam.groupby([bin_col, 'global_diameter_bin'], observed=False)
    .size()
    .unstack(fill_value=0)
)
global_diameter_distr.columns = [f'Diameter_{c}' for c in global_diameter_distr.columns]

stats = stats.join(global_area_distr, how='left')
stats = stats.join(global_diameter_distr, how='left')

summary = stats.agg(['mean', 'std', 'sem'])
summary.index = ['Mean', 'Std', 'SEM']
stats_with_summary = pd.concat([stats, summary])

# -------------------------
# Global diameter percentiles and fit
# -------------------------
global_d10, global_d50, global_d90 = compute_percentiles(df['Diameter'])

fit_diameter = df['Diameter']
fit_diameter = fit_diameter[fit_diameter > 0]

mu, sigma = lognormal_fit(fit_diameter)

x_fit = np.linspace(fit_diameter.min(), fit_diameter.max(), 300) if len(fit_diameter) > 0 else np.array([])
pdf_fit = lognormal_pdf(x_fit, mu, sigma) if len(x_fit) > 0 else np.array([])

fit_curve_df = pd.DataFrame({
    'Diameter': x_fit,
    'PDF': pdf_fit
})

# Global density from diameter
hist_counts, hist_edges = np.histogram(fit_diameter, bins=global_diameter_bins)
hist_centers = (hist_edges[:-1] + hist_edges[1:]) / 2
hist_widths = np.diff(hist_edges)
hist_density = hist_counts / (hist_counts.sum() * hist_widths) if hist_counts.sum() > 0 else np.zeros_like(hist_centers)

global_density_df = pd.DataFrame({
    'Bin Center': hist_centers,
    'Count': hist_counts,
    'Density': hist_density
})

# -------------------------
# Output file
# -------------------------
if args.o:
    output_filename = args.o if args.o.lower().endswith('.xlsx') else args.o + '.xlsx'
else:
    output_filename = f'{args.input[:-4]}_stats{n_bin}.xlsx'

# -------------------------
# Write workbook
# -------------------------
with pd.ExcelWriter(output_filename, engine='xlsxwriter') as writer:
    workbook = writer.book

    # Summary sheet
    stats_with_summary.to_excel(writer, sheet_name='Summary')
    summary_ws = writer.sheets['Summary']

    metrics_col = stats_with_summary.shape[1] + 3
    summary_ws.write(0, metrics_col, 'Global Diameter Metrics')
    summary_ws.write(1, metrics_col, 'D10')
    summary_ws.write(1, metrics_col + 1, global_d10)
    summary_ws.write(2, metrics_col, 'D50')
    summary_ws.write(2, metrics_col + 1, global_d50)
    summary_ws.write(3, metrics_col, 'D90')
    summary_ws.write(3, metrics_col + 1, global_d90)
    summary_ws.write(4, metrics_col, 'Lognormal mu')
    summary_ws.write(4, metrics_col + 1, mu)
    summary_ws.write(5, metrics_col, 'Lognormal sigma')
    summary_ws.write(5, metrics_col + 1, sigma)

    # Fit sheet
    fit_curve_df.to_excel(writer, sheet_name='Global_LogNormal_Diameter', index=False)
    fit_ws = writer.sheets['Global_LogNormal_Diameter']
    fit_ws.write(0, 3, 'mu')
    fit_ws.write(0, 4, mu)
    fit_ws.write(1, 3, 'sigma')
    fit_ws.write(1, 4, sigma)
    fit_ws.write(2, 3, 'D10')
    fit_ws.write(2, 4, global_d10)
    fit_ws.write(3, 3, 'D50')
    fit_ws.write(3, 4, global_d50)
    fit_ws.write(4, 3, 'D90')
    fit_ws.write(4, 4, global_d90)

    # Density sheet
    global_density_df.to_excel(writer, sheet_name='Global_Density_Diameter', index=False)
    gd_ws = writer.sheets['Global_Density_Diameter']

    chart_global = workbook.add_chart({'type': 'scatter', 'subtype': 'straight'})
    chart_global.add_series({
        'name': 'Global Density',
        'categories': ['Global_Density_Diameter', 1, 0, len(global_density_df), 0],
        'values':     ['Global_Density_Diameter', 1, 2, len(global_density_df), 2],
    })
    if len(fit_curve_df) > 0:
        chart_global.add_series({
            'name': 'Log-normal Fit',
            'categories': ['Global_LogNormal_Diameter', 1, 0, len(fit_curve_df), 0],
            'values':     ['Global_LogNormal_Diameter', 1, 1, len(fit_curve_df), 1],
        })
    chart_global.set_title({'name': 'Global Diameter Density and Log-normal Fit'})
    chart_global.set_x_axis({'name': 'Diameter'})
    chart_global.set_y_axis({'name': 'Density'})
    gd_ws.insert_chart(1, 5, chart_global)

    # Sector sheets
    total_area_all = df['Area'].sum()

    for sector in range(n_bin):
        sector_df = df[df[bin_col] == sector].copy()
        sheet_name = f'Sector_{sector}'

        sector_df.to_excel(writer, sheet_name=sheet_name, index=False)
        ws = writer.sheets[sheet_name]

        total_area = sector_df['Area'].sum()
        particle_count = len(sector_df)
        area_fraction = total_area / total_area_all if total_area_all > 0 else np.nan
        d10, d50, d90 = compute_percentiles(sector_df['Diameter'])
        mu_s, sigma_s = lognormal_fit(sector_df['Diameter'])

        info_col = len(sector_df.columns) + 2
        ws.write(0, info_col, 'Sector Summary')
        ws.write(1, info_col, 'Particle Count')
        ws.write(1, info_col + 1, particle_count)
        ws.write(2, info_col, 'Total Area')
        ws.write(2, info_col + 1, total_area)
        ws.write(3, info_col, 'Area Fraction')
        ws.write(3, info_col + 1, area_fraction)
        ws.write(4, info_col, 'D10 Diameter')
        ws.write(4, info_col + 1, d10)
        ws.write(5, info_col, 'D50 Diameter')
        ws.write(5, info_col + 1, d50)
        ws.write(6, info_col, 'D90 Diameter')
        ws.write(6, info_col + 1, d90)
        ws.write(7, info_col, 'Lognormal mu (Diameter)')
        ws.write(7, info_col + 1, mu_s)
        ws.write(8, info_col, 'Lognormal sigma (Diameter)')
        ws.write(8, info_col + 1, sigma_s)

        # -------------------------
        # Global Area bins
        # -------------------------
        start_row = len(sector_df) + 3
        sector_area_global = sector_df.copy()
        if args.log_bins:
            sector_area_global = sector_area_global[sector_area_global['Area'] > 0].copy()

        global_area_table = make_distribution_table(
            sector_area_global['Area'],
            global_area_bins,
            'Global Area Bin'
        )
        global_area_table.to_excel(writer, sheet_name=sheet_name, startrow=start_row, index=False)

        chart1 = workbook.add_chart({'type': 'column'})
        if len(global_area_table) > 0:
            chart1.add_series({
                'name': f'Sector {sector} Global Area',
                'categories': [sheet_name, start_row + 1, 0, start_row + len(global_area_table), 0],
                'values':     [sheet_name, start_row + 1, 1, start_row + len(global_area_table), 1],
            })
        chart1.set_title({'name': f'Sector {sector} Global Area Bins'})
        chart1.set_x_axis({'name': 'Area Bin'})
        chart1.set_y_axis({'name': 'Count'})
        ws.insert_chart(start_row, 3, chart1)

        # -------------------------
        # Local Area bins
        # -------------------------
        start_row2 = start_row + len(global_area_table) + 18
        local_area_bins = make_bins(
            sector_df['Area'],
            args.n_bin_size,
            use_log=args.log_bins
        )
        local_area_table = make_distribution_table(
            sector_df['Area'],
            local_area_bins,
            'Local Area Bin'
        )
        local_area_table.to_excel(writer, sheet_name=sheet_name, startrow=start_row2, index=False)

        chart2 = workbook.add_chart({'type': 'column'})
        if len(local_area_table) > 0:
            chart2.add_series({
                'name': f'Sector {sector} Local Area',
                'categories': [sheet_name, start_row2 + 1, 0, start_row2 + len(local_area_table), 0],
                'values':     [sheet_name, start_row2 + 1, 1, start_row2 + len(local_area_table), 1],
            })
        chart2.set_title({'name': f'Sector {sector} Local Area Bins'})
        chart2.set_x_axis({'name': 'Area Bin'})
        chart2.set_y_axis({'name': 'Count'})
        ws.insert_chart(start_row2, 3, chart2)

        # -------------------------
        # Global Diameter bins
        # -------------------------
        start_row3 = start_row2 + len(local_area_table) + 18
        sector_diam_global = sector_df.copy()
        if args.log_bins:
            sector_diam_global = sector_diam_global[sector_diam_global['Diameter'] > 0].copy()

        global_diameter_table = make_distribution_table(
            sector_diam_global['Diameter'],
            global_diameter_bins,
            'Global Diameter Bin'
        )
        global_diameter_table.to_excel(writer, sheet_name=sheet_name, startrow=start_row3, index=False)

        chart3 = workbook.add_chart({'type': 'column'})
        if len(global_diameter_table) > 0:
            chart3.add_series({
                'name': f'Sector {sector} Global Diameter',
                'categories': [sheet_name, start_row3 + 1, 0, start_row3 + len(global_diameter_table), 0],
                'values':     [sheet_name, start_row3 + 1, 1, start_row3 + len(global_diameter_table), 1],
            })
        chart3.set_title({'name': f'Sector {sector} Global Diameter Bins'})
        chart3.set_x_axis({'name': 'Diameter Bin'})
        chart3.set_y_axis({'name': 'Count'})
        ws.insert_chart(start_row3, 3, chart3)

        # -------------------------
        # Local Diameter bins
        # -------------------------
        start_row4 = start_row3 + len(global_diameter_table) + 18
        local_diameter_bins = make_bins(
            sector_df['Diameter'],
            args.n_bin_size,
            use_log=args.log_bins
        )
        local_diameter_table = make_distribution_table(
            sector_df['Diameter'],
            local_diameter_bins,
            'Local Diameter Bin'
        )
        local_diameter_table.to_excel(writer, sheet_name=sheet_name, startrow=start_row4, index=False)

        chart4 = workbook.add_chart({'type': 'column'})
        if len(local_diameter_table) > 0:
            chart4.add_series({
                'name': f'Sector {sector} Local Diameter',
                'categories': [sheet_name, start_row4 + 1, 0, start_row4 + len(local_diameter_table), 0],
                'values':     [sheet_name, start_row4 + 1, 1, start_row4 + len(local_diameter_table), 1],
            })
        chart4.set_title({'name': f'Sector {sector} Local Diameter Bins'})
        chart4.set_x_axis({'name': 'Diameter Bin'})
        chart4.set_y_axis({'name': 'Count'})
        ws.insert_chart(start_row4, 3, chart4)

print(f'\nResults saved to {output_filename}')
