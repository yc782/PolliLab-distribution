"""Whole-tablet ImageJ particle statistics and sector-based spatial CVs.

Input convention follows TabletDistributionFinal.py: first CSV data row stores
center coordinates (X,Y), followed by particle rows; last four lines are footer.
Areas and coordinates retain the units supplied by ImageJ.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def cv_percent(values):
    values = np.asarray(values, dtype=float)
    mean = values.mean()
    return float(100 * values.std(ddof=1) / mean) if len(values) > 1 and mean > 0 else np.nan


def make_edges(values, n_bins, logarithmic=False, minimum=None, maximum=None):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if logarithmic:
        values = values[values > 0]
    if not len(values):
        return None
    lo = float(values.min()) if minimum is None else minimum
    hi = float(values.max()) if maximum is None else maximum
    if logarithmic and lo <= 0:
        raise ValueError('Logarithmic bin minimum must be positive')
    if hi < lo:
        raise ValueError('Bin maximum must be >= minimum')
    if hi == lo:
        hi = lo * (1 + 1e-6) if logarithmic else lo + max(abs(lo) * 1e-6, 1e-9)
    return np.geomspace(lo, hi, n_bins + 1) if logarithmic else np.linspace(lo, hi, n_bins + 1)


def distribution(values, edges):
    if edges is None:
        return pd.DataFrame(columns=['Bin lower', 'Bin upper', 'Count', 'Percent of particles'])
    counts, _ = np.histogram(values, bins=edges)
    total = len(values)
    return pd.DataFrame({
        'Bin lower': edges[:-1], 'Bin upper': edges[1:],
        'Count': counts, 'Percent of particles': 100 * counts / total if total else 0,
    })


def main():
    parser = argparse.ArgumentParser(description='Whole-tablet particle metrics and spatial homogeneity')
    parser.add_argument('input', help='ImageJ particle CSV (same format as original script)')
    parser.add_argument('-o', help='Output XLSX filename')
    parser.add_argument('--n_bin', type=int, default=16, help='Number of equal-angle sectors (default 16)')
    parser.add_argument('--n_bin_size', type=int, default=10, help='Number of area/diameter histogram bins')
    parser.add_argument('--min_area', type=float, default=None)
    parser.add_argument('--max_area', type=float, default=None)
    parser.add_argument('--min_diameter', type=float, default=None)
    parser.add_argument('--max_diameter', type=float, default=None)
    parser.add_argument('--log_bins', action='store_true')
    args = parser.parse_args()
    if args.n_bin < 2 or args.n_bin_size < 1:
        parser.error('--n_bin must be >= 2 and --n_bin_size >= 1')
    if Path(args.input).suffix.lower() != '.csv':
        parser.error('Input must be a CSV file')

    raw = pd.read_csv(args.input, skipfooter=4, engine='python')
    if not {'X', 'Y', 'Area'}.issubset(raw.columns) or len(raw) < 2:
        raise ValueError('CSV must have X, Y, Area columns, a center row, and particle rows')
    center_x = float(raw.iloc[0]['X'])
    center_y = float(raw.iloc[0]['Y'])
    df = raw.iloc[1:].copy()
    for col in ('X', 'Y', 'Area'):
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['X', 'Y', 'Area'])
    df = df[df['Area'] > 0].copy()
    if df.empty:
        raise ValueError('No valid positive-area particles found')

    df['Diameter'] = np.sqrt(4 * df['Area'] / np.pi)
    theta = np.arctan2(df['Y'] - center_y, df['X'] - center_x)
    df['Sector'] = (np.floor((theta + np.pi) / (2 * np.pi) * args.n_bin).astype(int) % args.n_bin)

    # Reindex ensures sectors containing zero particles contribute zero to spatial CV.
    sectors = pd.RangeIndex(args.n_bin)
    sector_counts = df.groupby('Sector').size().reindex(sectors, fill_value=0)
    sector_areas = df.groupby('Sector')['Area'].sum().reindex(sectors, fill_value=0.0)
    area = df['Area'].to_numpy()
    diameter = df['Diameter'].to_numpy()
    log_diameter = np.log(diameter)
    mu = float(log_diameter.mean())
    sigma = float(log_diameter.std(ddof=1)) if len(log_diameter) > 1 else np.nan

    metrics = pd.DataFrame([
        ('Number of sectors', args.n_bin),
        ('Total particle count', len(df)),
        ('Particle count CV across sectors (%)', cv_percent(sector_counts)),
        ('Particle area-sum CV across sectors (%)', cv_percent(sector_areas)),
        ('Mean particle area', float(area.mean())),
        ('Median particle area', float(np.median(area))),
        ('Particle area SD', float(area.std(ddof=1)) if len(area) > 1 else np.nan),
        ('Mean equivalent-circle particle diameter', float(diameter.mean())),
        ('Equivalent-circle particle diameter SD', float(diameter.std(ddof=1)) if len(diameter) > 1 else np.nan),
        ('Lognormal diameter mu', mu),
        ('Lognormal diameter sigma', sigma),
    ], columns=['Metric', 'Value'])

    area_edges = make_edges(area, args.n_bin_size, args.log_bins, args.min_area, args.max_area)
    diameter_edges = make_edges(diameter, args.n_bin_size, args.log_bins, args.min_diameter, args.max_diameter)
    area_dist = distribution(area, area_edges)
    diameter_dist = distribution(diameter, diameter_edges)
    fit_x = np.linspace(diameter.min(), diameter.max(), 300)
    if np.isfinite(sigma) and sigma > 0:
        fit_y = np.exp(-0.5 * ((np.log(fit_x) - mu) / sigma) ** 2) / (fit_x * sigma * np.sqrt(2 * np.pi))
    else:
        fit_y = np.full_like(fit_x, np.nan)
    fit_df = pd.DataFrame({'Diameter': fit_x, 'Lognormal PDF': fit_y})

    notes = pd.DataFrame({'Methodological notes': [
        'Input is ImageJ particle CSV; first data row is tablet center and last four CSV lines are skipped, matching original script.',
        'Sector assignments use particle centroid (X,Y) and equal-angle radial sectors.',
        'Count CV and area-sum CV use sample standard deviation (ddof=1), including empty sectors as zeros.',
        'Area-sum CV is a proxy for coverage CV only if sectors have equal valid tablet area and particles do not substantially cross boundaries.',
        'Particles spanning boundaries are assigned entirely to their centroid sector; pixel-level coverage CV requires the segmented image.',
        'ImageJ calibration determines area and diameter units; equivalent-circle diameter is calculated from area.',
        'Mean diameter and diameter SD are calculated across all retained particles; SD uses ddof=1 (sample SD), not the lognormal fit sigma.',
        'Histogram percentages use all detected particles as denominator; values outside custom bin limits are not counted in bins.',
        'Spatial CV measures between-sector variation, not within-sector clustering or independent-batch variability.',
    ]})

    output = Path(args.o) if args.o else Path(args.input).with_name(Path(args.input).stem + '_summary.xlsx')
    if output.suffix.lower() != '.xlsx':
        output = output.with_suffix('.xlsx')
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        metrics.to_excel(writer, sheet_name='Summary', index=False)
        area_dist.to_excel(writer, sheet_name='Particle_Area_Distribution', index=False)
        diameter_dist.to_excel(writer, sheet_name='Diameter_Distribution', index=False)
        fit_df.to_excel(writer, sheet_name='Diameter_Lognormal_Fit', index=False)
        notes.to_excel(writer, sheet_name='Methods_Notes', index=False)
        book = writer.book
        for sheet, table, x_label in [
            ('Particle_Area_Distribution', area_dist, 'Particle area bin'),
            ('Diameter_Distribution', diameter_dist, 'Diameter bin'),
        ]:
            ws = writer.sheets[sheet]
            ws.set_column('A:B', 18)
            ws.set_column('C:D', 22)
            chart = book.add_chart({'type': 'column'})
            chart.add_series({
                'name': 'Particle count',
                'categories': [sheet, 1, 0, len(table), 0],
                'values': [sheet, 1, 2, len(table), 2],
            })
            chart.set_title({'name': x_label + ' distribution'})
            chart.set_x_axis({'name': x_label + ' (lower edge)'})
            chart.set_y_axis({'name': 'Particle count'})
            ws.insert_chart('F2', chart)
        writer.sheets['Summary'].set_column('A:A', 44)
        writer.sheets['Summary'].set_column('B:B', 22)
        writer.sheets['Methods_Notes'].set_column('A:A', 115)
    print(f'Results saved to {output}')


if __name__ == '__main__':
    main()
