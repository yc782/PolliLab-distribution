"""Sensitivity analysis for centroid-based QCL-IR particle spatial CVs.

Input: same ImageJ CSV convention as TabletDistributionUpdated.py:
first data row = tablet center X,Y; last four lines = footer (by default).
All particles are assigned to sectors by their ImageJ centroid.

Example:
  python TabletSensitivityAnalysis.py tablet.csv
  python TabletSensitivityAnalysis.py tablet.csv --sectors 8 16 32 --rotations 0 0.25 0.5 0.75 -o sensitivity.xlsx
  python TabletSensitivityAnalysis.py tablet.csv --min-areas 0 10 25

Rotations are FRACTIONS of one sector width, not degrees; 0.5 means half a sector.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def cv_percent(values):
    values = np.asarray(values, dtype=float)
    if values.size < 2 or not np.isfinite(values).all():
        return np.nan
    mean = values.mean()
    return float(100 * values.std(ddof=1) / mean) if mean > 0 else np.nan


def load_particles(path, footer_lines):
    raw = pd.read_csv(path, skipfooter=footer_lines, engine='python')
    if not {'X', 'Y', 'Area'}.issubset(raw.columns) or len(raw) < 2:
        raise ValueError('Expected X, Y, Area columns, a first data row containing tablet center, and particle rows.')
    center_x = float(raw.iloc[0]['X'])
    center_y = float(raw.iloc[0]['Y'])
    df = raw.iloc[1:].copy()
    for col in ['X', 'Y', 'Area']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['X', 'Y', 'Area'])
    df = df[df['Area'] > 0].copy()
    if df.empty:
        raise ValueError('No positive-area particles were found.')
    df['theta'] = np.arctan2(df['Y'] - center_y, df['X'] - center_x)
    return df


def analyze(df, n_sectors, rotation_fraction, minimum_area):
    selected = df[df['Area'] >= minimum_area]
    angles = selected['theta'].to_numpy()
    areas = selected['Area'].to_numpy()
    # Rotating boundaries by a fraction of sector width; fixed origin at -pi.
    sector = np.floor((angles + np.pi) / (2 * np.pi) * n_sectors - rotation_fraction).astype(int) % n_sectors
    counts = np.bincount(sector, minlength=n_sectors)
    summed_area = np.bincount(sector, weights=areas, minlength=n_sectors)
    return {
        'Sectors': n_sectors,
        'Rotation (sector fractions)': rotation_fraction,
        'Rotation (degrees)': rotation_fraction * 360.0 / n_sectors,
        'Minimum particle area': minimum_area,
        'Particles retained': len(selected),
        'Total retained particle area': float(areas.sum()),
        'Particle count CV (%)': cv_percent(counts),
        'Particle area-sum CV (%)': cv_percent(summed_area),
        'Empty sectors': int((counts == 0).sum()),
    }


def main():
    p = argparse.ArgumentParser(description='Sensitivity analysis of QCL-IR sector-wise particle count and area-sum CVs')
    p.add_argument('input', help='ImageJ CSV (same format as TabletDistributionUpdated.py)')
    p.add_argument('-o', '--output', help='Output Excel .xlsx filename')
    p.add_argument('--sectors', type=int, nargs='+', default=[8, 16, 32], help='Numbers of radial sectors (default: 8 16 32)')
    p.add_argument('--rotations', type=float, nargs='+', default=[0, 0.25, 0.5, 0.75],
                   help='Boundary rotations as fractions of ONE sector width (default: 0 0.25 0.5 0.75)')
    p.add_argument('--min-areas', type=float, nargs='+', default=[0],
                   help='Minimum particle-area inclusion thresholds, in ImageJ area units (default: 0)')
    p.add_argument('--footer-lines', type=int, default=4, help='CSV footer lines to skip (default: 4)')
    args = p.parse_args()
    if any(n < 2 for n in args.sectors):
        p.error('All sector counts must be at least 2.')
    if any(not 0 <= r < 1 for r in args.rotations):
        p.error('All rotation fractions must be in [0, 1).')
    if any(x < 0 for x in args.min_areas):
        p.error('Minimum particle areas must be nonnegative.')
    if args.footer_lines < 0:
        p.error('--footer-lines must be nonnegative.')
    df = load_particles(args.input, args.footer_lines)
    rows = [analyze(df, n, r, m)
            for m in sorted(set(args.min_areas))
            for n in sorted(set(args.sectors))
            for r in sorted(set(args.rotations))]
    detail = pd.DataFrame(rows)
    grouping = ['Minimum particle area', 'Sectors']
    grouped = detail.groupby(grouping, as_index=False)
    summary = grouped.agg(
        Rotations_tested=('Rotation (degrees)', 'count'),
        Count_CV_mean=('Particle count CV (%)', 'mean'),
        Count_CV_min=('Particle count CV (%)', 'min'),
        Count_CV_max=('Particle count CV (%)', 'max'),
        Count_CV_SD=('Particle count CV (%)', 'std'),
        Area_CV_mean=('Particle area-sum CV (%)', 'mean'),
        Area_CV_min=('Particle area-sum CV (%)', 'min'),
        Area_CV_max=('Particle area-sum CV (%)', 'max'),
        Area_CV_SD=('Particle area-sum CV (%)', 'std'),
        Empty_sectors_max=('Empty sectors', 'max'),
    )
    # Baseline uses 16 sectors and unrotated boundaries, matching existing script.
    baseline = detail[(detail['Sectors'] == 16) & (np.isclose(detail['Rotation (sector fractions)'], 0)) &
                      (np.isclose(detail['Minimum particle area'], 0))]
    baseline_note = baseline[['Particle count CV (%)', 'Particle area-sum CV (%)']].copy()
    baseline_note.insert(0, 'Baseline definition', '16 sectors; no rotation; all particles')
    if baseline_note.empty:
        baseline_note = pd.DataFrame({'Baseline definition': ['Not included in requested parameter grid']})
    notes = pd.DataFrame({'Interpretation and limitations': [
        'Rotations are fractions of the width of ONE sector; 0.25 is one-quarter of that sector angle.',
        'Area CV is CV of SUMMED ImageJ particle area per sector, not CV of individual particle areas.',
        'Both CVs use sample SD (ddof=1) to match TabletDistributionUpdated.py.',
        'Empty sectors contribute zero to both CVs; high sector counts may increase CV due to sparse particle counts.',
        'The same image is repeatedly reanalyzed; settings are not independent biological or batch replicates.',
        'Particles crossing sector boundaries are wholly assigned by centroid; this is not pixel-based coverage.',
        'Area thresholds change the analyzed particle population and total retained area; use only scientifically justified thresholds.',
        'Sector-count effects and boundary-rotation effects are distinct: compare rotations WITHIN each sector count.',
        'A robust process ranking should ideally persist across reasonable sector counts and boundary rotations.',
        'For exact red-pixel coverage CV and segmentation-threshold sensitivity, analyze original segmented images.',
        'Sector CV describes between-sector spatial variation, not clustering within sectors.',
    ]})
    output = Path(args.output) if args.output else Path(args.input).with_name(Path(args.input).stem + '_sensitivity.xlsx')
    if output.suffix.lower() != '.xlsx':
        output = output.with_suffix('.xlsx')
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        detail.to_excel(writer, sheet_name='All_Settings', index=False)
        summary.to_excel(writer, sheet_name='Rotation_Summary', index=False)
        baseline_note.to_excel(writer, sheet_name='Baseline', index=False)
        notes.to_excel(writer, sheet_name='Methods_Notes', index=False)
        for name in ['All_Settings', 'Rotation_Summary', 'Baseline']:
            ws = writer.sheets[name]
            ws.freeze_panes(1, 0)
            ws.set_column(0, 15, 23)
        writer.sheets['Methods_Notes'].set_column(0, 0, 115)
    print(f'Sensitivity analysis saved to {output}')
    print(f'Analyzed {len(df)} particles across {len(detail)} settings.')


if __name__ == '__main__':
    main()
