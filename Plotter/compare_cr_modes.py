import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt

import argparse
import re
import math
import sys, os
from typing import List, Dict, Tuple
import numpy as np
from matplotlib.colors import to_rgba, to_rgb

# ====================================================================
# YODA PARSING AND CALCULATION FUNCTIONS (unchanged)
# ====================================================================

def read_file_content(filename: str) -> str:
    """Reads content from a specified file, handling FileNotFoundError."""
    try:
        with open(filename, "r") as f:
            return f.read()
    except FileNotFoundError:
        print(f"Error: Component file '{filename}' not found.")
        return ""

def parse_yoda_values_and_scale(yoda_data: str, path: str) -> Tuple[List[float], float]:
    """Parse a YODA object at `path` and return a per-bin array plus an overall scale."""
    escaped_path = re.escape(path)

    # Try parsing YODA_ESTIMATE1D_V3 blocks
    est_pattern = rf'BEGIN YODA_ESTIMATE1D_V3 {escaped_path}.*?ScaledBy: ([\d\.\+\-e]+).*?Edges\(A1\): \[(.*?)\].*?# value.*?\n(.*?)\nEND YODA_ESTIMATE1D_V3'
    est_match = re.search(est_pattern, yoda_data, re.DOTALL)

    if est_match:
        scale_factor = float(est_match.group(1))
        edges = [float(e) for e in est_match.group(2).split(',')]
        data_block = est_match.group(3).strip()

        values = []
        for line in data_block.split('\n'):
            parts = line.split()
            if parts[0].lower() == 'nan':
                continue
            try:
                values.append(float(parts[0]) * scale_factor)
            except ValueError:
                continue

        return values, scale_factor

    # Try parsing YODA_HISTO1D_V3 blocks
    histo_pattern = rf'BEGIN YODA_HISTO1D_V3 {escaped_path}.*?Edges\(A1\): \[(.*?)\].*?# sumW.*?\n(.*?)\nEND YODA_HISTO1D_V3'
    histo_match = re.search(histo_pattern, yoda_data, re.DOTALL)

    if histo_match:
        edges = [float(e) for e in histo_match.group(1).split(',')]
        data_block = histo_match.group(2).strip()

        values = []
        for line in data_block.split('\n'):
            parts = line.split()
            if len(parts) >= 1:
                try:
                    values.append(float(parts[0]))  # Use the first column (sumW)
                except ValueError:
                    continue

        return values, 1.0  # No explicit scale in Histo1D blocks

    return [], 1.0

def get_yoda_template_info(yoda_data: str, path: str) -> Tuple[List[float], List[str], str]:
    """Extracts the bin edges, error labels, and title from a YODA object for template."""
    escaped_path = re.escape(path)

    # Try parsing YODA_ESTIMATE1D_V3 blocks
    est_pattern = rf'BEGIN YODA_ESTIMATE1D_V3 {escaped_path}.*?Title:\s*(.*?)\n.*?Edges\(A1\): \[(.*?)\].*?ErrorLabels: \[(.*?)\]'
    est_match = re.search(est_pattern, yoda_data, re.DOTALL)

    if est_match:
        title = est_match.group(1).strip()
        edges = [float(e) for e in est_match.group(2).split(',')]
        error_labels = [label.strip('"') for label in est_match.group(3).split(',')]
        return edges, error_labels, title

    # Try parsing YODA_HISTO1D_V3 blocks
    histo_pattern = rf'BEGIN YODA_HISTO1D_V3 {escaped_path}.*?Title:\s*(.*?)\n.*?Edges\(A1\): \[(.*?)\]'
    histo_match = re.search(histo_pattern, yoda_data, re.DOTALL)

    if histo_match:
        title = histo_match.group(1).strip()
        edges = [float(e) for e in histo_match.group(2).split(',')]
        return edges, [], title

    # Default return if no match is found
    return [], [], ""

def calculate_normalized_data(
    Pwave_data: List[float], 
    Swave_data: List[float], 
    Singlet_data: List[float], 
    Octet_data: List[float]
) -> Tuple[Dict[str, List[float]], float]:
    """Calculates the normalized data for all six quantities, with per-channel normalization."""
    component_count = len(Pwave_data)
    if component_count == 0:
        print("Error: No data provided for normalization.")
        return {}, 0.0

    # Per-channel normalization factors
    norm_Pwave = sum(Pwave_data) if sum(Pwave_data) != 0 else 1.0
    norm_Swave = sum(Swave_data) if sum(Swave_data) != 0 else 1.0
    norm_Singlet = sum(Singlet_data) if sum(Singlet_data) != 0 else 1.0
    norm_Octet = sum(Octet_data) if sum(Octet_data) != 0 else 1.0
    norm_FeedDown = sum([Pwave_data[i] + Swave_data[i] for i in range(component_count)])
    if norm_FeedDown == 0:
        norm_FeedDown = 1.0
    Total_Unnormalized_Bins = [
        Pwave_data[i] + Swave_data[i] + Singlet_data[i] + Octet_data[i]
        for i in range(component_count)
    ]
    norm_Total = sum(Total_Unnormalized_Bins) if sum(Total_Unnormalized_Bins) != 0 else 1.0

    Normalized_Data_Dict = {
        "Pwave": [val / norm_Pwave for val in Pwave_data],
        "Swave": [val / norm_Swave for val in Swave_data],
        "Singlet": [val / norm_Singlet for val in Singlet_data],
        "Octet": [val / norm_Octet for val in Octet_data],
        "FeedDown": [],
        "Total": []
    }

    for i in range(component_count):
        feeddown_unnorm_bin = Pwave_data[i] + Swave_data[i]
        Normalized_Data_Dict["FeedDown"].append(feeddown_unnorm_bin / norm_FeedDown)
        Normalized_Data_Dict["Total"].append(Total_Unnormalized_Bins[i] / norm_Total)

    print("Debugging normalization factors:")
    print(f"  Pwave: {norm_Pwave}, Swave: {norm_Swave}, Singlet: {norm_Singlet}, Octet: {norm_Octet}, FeedDown: {norm_FeedDown}, Total: {norm_Total}")
    print("Debugging Total Values:")
    for i, value in enumerate(Normalized_Data_Dict["Total"]):
        print(f"Bin {i}: {value}")

    # Store unnormalized bins and normalization factors for error propagation
    Unnorm_Bins = {
        "Pwave": Pwave_data,
        "Swave": Swave_data,
        "Singlet": Singlet_data,
        "Octet": Octet_data,
        "FeedDown": [Pwave_data[i] + Swave_data[i] for i in range(component_count)],
        "Total": Total_Unnormalized_Bins
    }
    Norm_Factors = {
        "Pwave": norm_Pwave,
        "Swave": norm_Swave,
        "Singlet": norm_Singlet,
        "Octet": norm_Octet,
        "FeedDown": norm_FeedDown,
        "Total": norm_Total
    }
    return Normalized_Data_Dict, Norm_Factors, Unnorm_Bins

def create_yoda_output_block(path: str, data: List[float], edges: List[float], title: str = "") -> str:
    """
    Formats calculated data into the requested YODA_HISTO1D_V3 block structure
    using placeholders for unknown statistical moments.
    """
    if not data or len(edges) != len(data) + 1:
        return ""

    sumw_total = sum(data)
    edges_str = ', '.join(f"{e:.6e}" for e in edges)

    output = f"BEGIN YODA_HISTO1D_V3 {path}\n"
    output += f"Path: {path}\n"
    output += "ScaledBy: 1.00000000000000000e+00\n"
    output += f"Title: {title}\n"
    output += "Type: Histo1D\n"
    output += "---\n"
    output += f"# Mean: {0.0:.6e}\n"
    output += f"# Integral: {sumw_total:.6e}\n"
    output += f"Edges(A1): [{edges_str}]\n"
    output += "# ID\t ID\t sumw\t sumw2\t sumwx\t sumwx2\t numEntries\n"
    output += f"Total\tTotal\t{sumw_total:.6e}\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\t{len(data):.1f}e+00\n"
    output += f"Underflow\tUnderflow\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\n"
    output += f"Overflow\tOverflow\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\n"
    output += "# xlow\t xhigh\t sumw\t sumw2\t sumwx\t sumwx2\t numEntries\n"

    for i in range(len(data)):
        xlow = edges[i]
        xhigh = edges[i + 1]
        sumw = data[i]
        if math.isnan(sumw):
            continue
        output += (
            f"{xlow:.6e}\t{xhigh:.6e}\t{sumw:.6e}\t"
            f"{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\t{0.0:.6e}\n"
        )

    output += "END YODA_HISTO1D_V3\n"
    return output


def write_normalized_yoda_outputs(
    edges: List[float],
    cr_mode_data: Dict[str, Dict[str, List[float]]],
    output_dir: str,
    system_value: str,
    energy_value: str,
    ldme_value: str,
    ptHat_value: str,
    event_value: str,
    rivet_config_value: str,
) -> None:
    """Write normalized channel outputs as YODA_HISTO1D files for each CR mode."""
    if not edges or not cr_mode_data:
        print("Skipping normalized YODA output: No valid data found.")
        return

    os.makedirs(output_dir, exist_ok=True)
    channels = ["Pwave", "Swave", "Singlet", "Octet", "FeedDown", "Total"]
    output_path_base = "/MyAnalysis/Normalized"

    for cr_mode in sorted(cr_mode_data.keys()):
        normed = cr_mode_data[cr_mode].get("normalized", {})
        for channel in channels:
            values = normed.get(channel, [])
            if not values:
                continue

            yoda_path = f"{output_path_base}/{channel}/CR{cr_mode}"
            yoda_title = f"{channel} Normalized CR{cr_mode}"
            yoda_block = create_yoda_output_block(yoda_path, values, edges, yoda_title)
            if not yoda_block:
                continue

            filename = (
                f"normalized_{channel.lower()}_{system_value}_{energy_value}_"
                f"ptHat{ptHat_value}_LDMEfac{ldme_value}_CR{cr_mode}_{event_value}_"
                f"{rivet_config_value}.yoda"
            )
            file_path = os.path.join(output_dir, filename)
            with open(file_path, "w") as f:
                f.write(yoda_block)

            print(f"Normalized YODA saved: {file_path}")


def resolve_hist_path(rivet_config_value: str) -> str:
    """Return the raw histogram path associated with the selected Rivet config."""
    hist_paths = {
        "LHCb_default": "/RAW/LHCB_2017_I1509507/d01-x01-y01",
        "LHCb_CMSCut_R3": "/RAW/LHCB_2017_I1509507_CMSCut_R3/frag_prompt",
        "LHCb_CMSCut_R4": "/RAW/LHCB_2017_I1509507_CMSCut_R4/frag_prompt",
        "LHCb_FinerBin": "/RAW/LHCB_2017_I1509507_FinerBin/frag_prompt",
        "CMS_default": "/RAW/CMS_2022_I1870319/d01-x01-y01",
        "CMS_FineBin_R3": "/RAW/CMS_2022_I1870319_FineBin_R3/frag_fixed_bins",
        "CMS_FineBin_R4": "/RAW/CMS_2022_I1870319_FineBin_R4/frag_fixed_bins",
        "CMS_20Bin_R3": "/RAW/CMS_2022_I1870319_20Bin_R3/frag_fixed_bins",
        "CMS_20Bin_R4": "/RAW/CMS_2022_I1870319_20Bin_R4/frag_fixed_bins",
        "CMS_50Bin_R3": "/RAW/CMS_2022_I1870319_50Bin_R3/frag_fixed_bins",
        "CMS_50Bin_R4": "/RAW/CMS_2022_I1870319_50Bin_R4/frag_fixed_bins",
        "CMS_100Bin_R3": "/RAW/CMS_2022_I1870319_100Bin_R3/frag_fixed_bins",
        "CMS_100Bin_R4": "/RAW/CMS_2022_I1870319_100Bin_R4/frag_fixed_bins",
        "CMS_LHCbcut_FineBin": "/RAW/CMS_2022_I1870319_LHCbcut_FineBin/frag_fixed_bins"
        
    }
    return hist_paths.get(rivet_config_value, hist_paths["LHCb_default"])


def resolve_x_limits(rivet_config_value: str) -> tuple[float, float]:
    """Return x-axis limits for the selected Rivet config."""
    if rivet_config_value.startswith("CMS_"):
        return 0.2, 1.0
    return 0.0, 1.0

def plot_comparison_for_cr_modes(
    edges: List[float],
    cr_mode_data: Dict[str, Dict[str, List[float]]],
    output_folder: str,
    system_value: str, Energy_value: str, ldme_value: str, ptHat_value: str, event_value: str, rivet_config_value: str):
    """Plots comparison for each channel across different CRModes."""
    # Debug print: Show 'Total' data for each CRMode before plotting
    print("\n[DEBUG] Data for 'Total' in each CRMode before plotting:")
    for cr_mode, data_dict in cr_mode_data.items():
        total_data = data_dict.get('Total', None)
        print("  CRMode {}: {}".format(cr_mode, total_data))
    """Plots comparison for each channel across different CRModes."""
    if not edges or not cr_mode_data:
        print("Skipping plot: No valid data found for plotting.")
        return

    centers = np.array([(edges[i] + edges[i+1]) / 2 for i in range(len(edges) - 1)])
    point_output_folder = os.path.join(output_folder, "point")
    hist_output_folder = os.path.join(output_folder, "hist")
    os.makedirs(point_output_folder, exist_ok=True)
    os.makedirs(hist_output_folder, exist_ok=True)

    plt.figure(figsize=(12, 8))

    styles = {
        "Octet":      {'color': 'red',     'label': 'Octet',                 'linewidth': 1.5},
        "Singlet":    {'color': 'blue',    'label': 'Singlet',               'linewidth': 1.5},
        "Pwave":      {'color': 'orange',  'label': 'Pwave',                 'linewidth': 1.5},
        "Swave":      {'color': 'gray',    'label': 'Swave',                 'linewidth': 1.5},
        "FeedDown":   {'color': 'green',   'label': 'FeedDown (P+S)',        'linewidth': 2.0},
        "Total":      {'color': 'black',   'label': 'All contribution',      'linewidth': 2.5}
    }
    channels = ["Pwave", "Swave", "Singlet", "Octet", "FeedDown", "Total"]

    for key in styles.keys():
        plt.figure(figsize=(8, 6))
        lines = []
        for cr_mode, cr_dict in cr_mode_data.items():
            normed = cr_dict["normalized"]
            norm_factors = cr_dict["norm_factors"]
            unnorm_bins = cr_dict["unnorm_bins"]
            if key in normed and normed[key]:
                base_color = 'red' if cr_mode == '0' else 'blue'
                color = base_color
                style = {k: v for k, v in styles[key].items() if k not in ['label', 'color']}
                # Error bars: sqrt(unnormalized bin) / normalization factor
                y = np.array(normed[key])
                yerr = np.array([np.sqrt(max(ub, 0)) / norm_factors[key] for ub in unnorm_bins[key]])
                line = plt.errorbar(
                    centers, y, yerr=yerr,
                    fmt='o',
                    color=color,
                    ecolor=color,
                    label=f"CRMode {cr_mode}",
                    markersize=7,
                    capsize=4,
                    elinewidth=3,
                    markeredgewidth=2,
                    linestyle='none',
                    **style
                )
                lines.append(line)

        plt.xlabel(r'z(J/$\psi$)', fontsize=12)
        plt.ylabel('Normalized entries', fontsize=12)
        plt.title(f'{system_value}, {Energy_value}, ptHat {ptHat_value}, LDMEfac {ldme_value}, {event_value}, {rivet_config_value}\n{key} Comparison', fontsize=13)
        if lines:
            plt.legend(loc='upper right', fontsize=10)
        # Annotate yields per CR mode beside legend (only for the plotted channel)
        ax = plt.gca()
        text_lines = []
        for cr_mode in sorted(cr_mode_data.keys()):
            vals = cr_mode_data[cr_mode]['unnorm_bins'].get(key, [])
            s = f"{float(np.sum(vals)):.1f}" if len(vals) > 0 else "0.0"
            text_lines.append(f"CR{cr_mode}: {s}")
        text_block = "\n".join(text_lines).rstrip()
        # place text left of the legend, right-aligned
        ax.text(0.78, 0.98, text_block, transform=ax.transAxes, fontsize=8, va='top', ha='right', family='monospace', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        plt.xlim(*resolve_x_limits(rivet_config_value))
        if rivet_config_value == "LHCb_default":
            plt.ylim(0, 0.4)
        plt.tight_layout()

        # Save the plot to the point-style folder.
        output_filename = f"{point_output_folder}/{key}_comparison.png"
        plt.savefig(output_filename, transparent=True)
        plt.close()
        print(f"Plot saved: {output_filename}")

    plt.figure(figsize=(8, 6))
    lines = []
    for cr_mode, cr_dict in cr_mode_data.items():
        normed = cr_dict["normalized"]
        norm_factors = cr_dict["norm_factors"]
        unnorm_bins = cr_dict["unnorm_bins"]
        total_data = normed.get('Total', None)
        if total_data:
            base_color = 'red' if cr_mode == '0' else 'blue'
            color = base_color
            y = np.array(total_data)
            yerr = np.array([np.sqrt(max(ub, 0)) / norm_factors['Total'] for ub in unnorm_bins['Total']])
            line = plt.errorbar(
                centers, y, yerr=yerr,
                fmt='o',
                color=color,
                ecolor=color,
                label="CRMode {}".format(cr_mode),
                markersize=10,
                capsize=4,
                elinewidth=3,
                markeredgewidth=2,
                linestyle='none',
                linewidth=2.5
            )
            lines.append(line)
    plt.xlabel(r'z(J/$\psi$)', fontsize=12)
    plt.ylabel('Normalized entries', fontsize=12)
    plt.title(f'{system_value}, {Energy_value}, ptHat {ptHat_value}, LDMEfac {ldme_value}, {event_value}, {rivet_config_value}\nAll Channels Comparison', fontsize=13)
    if lines:
        plt.legend(loc='upper right', fontsize=10)
    # Dynamic y-axis for Total plot: 1.4 * max(total + err) across CR modes
    max_y = 0.0
    for cr_mode, cr_dict in cr_mode_data.items():
        normed_local = cr_dict.get('normalized', {})
        if 'Total' in normed_local and normed_local['Total']:
            vals = np.array(normed_local['Total'])
            unnorm = cr_dict.get('unnorm_bins', {}).get('Total', [0.0] * len(vals))
            normf = cr_dict.get('norm_factors', {}).get('Total', 1.0) or 1.0
            yerr_local = np.array([np.sqrt(max(ub, 0)) / normf for ub in unnorm]) if len(unnorm) == len(vals) else np.zeros_like(vals)
            candidate = float(np.max(vals + yerr_local)) if vals.size else 0.0
            if candidate > max_y:
                max_y = candidate
    if max_y <= 0.0:
        max_y = 1.0
    plt.ylim(0, max_y * 1.4)
    plt.tight_layout()
    # Annotate yields next to legend for total plot (show all channels)
    ax = plt.gca()
    text_lines = []
    for cr_mode in sorted(cr_mode_data.keys()):
        text_lines.append(f"CR{cr_mode}:")
        for ch in channels:
            vals = cr_mode_data[cr_mode]['unnorm_bins'].get(ch, [])
            s = f"{float(np.sum(vals)):.1f}" if len(vals) > 0 else "0.0"
            text_lines.append(f" {ch[:3]}: {s}")
        text_lines.append("")
    text_block = "\n".join(text_lines).rstrip()
    # place text left of the legend, right-aligned
    ax.text(0.78, 0.98, text_block, transform=ax.transAxes, fontsize=8, va='top', ha='right', family='monospace', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

    plt.xlim(*resolve_x_limits(rivet_config_value))
    if rivet_config_value == "LHCb_default":
        plt.ylim(0, 0.4)

    total_output_filename = f"{point_output_folder}/comparison_CRModes.png"
    plt.savefig(total_output_filename, transparent=True)
    # plt.show()  # Suppressed for non-interactive backend
    print("\nComparison plot saved to '{}'".format(total_output_filename))

    # Histogram-style version of the same comparison set.
    for key in styles.keys():
        fig, ax = plt.subplots(figsize=(8, 6))
        plotted = False
        for cr_mode, cr_dict in cr_mode_data.items():
            normed = cr_dict["normalized"]
            if key in normed and normed[key]:
                color = 'red' if cr_mode == '0' else 'blue'
                ax.stairs(
                    np.array(normed[key]),
                    np.array(edges),
                    color=color,
                    label=f"CRMode {cr_mode}",
                    linewidth=2.0,
                    fill=False,
                )
                plotted = True

        ax.set_xlabel(r'z(J/$\psi$)', fontsize=12)
        ax.set_ylabel('Normalized entries', fontsize=12)
        ax.set_title(f'{system_value}, {Energy_value}, ptHat {ptHat_value}, LDMEfac {ldme_value}, {event_value}, {rivet_config_value}\n{key} Comparison (Histogram)', fontsize=13)
        if plotted:
            ax.legend(loc='upper right', fontsize=10)
        # Dynamic y-axis for histogram: 1.4 * max(normed values) for this channel
        max_y = 0.0
        for cr_mode, cr_dict in cr_mode_data.items():
            normed_local = cr_dict.get('normalized', {})
            if key in normed_local and normed_local[key]:
                vals = np.array(normed_local[key])
                if vals.size:
                    max_y = max(max_y, float(np.max(vals)))
        if max_y <= 0.0:
            max_y = 1.0
        ax.set_ylim(0, max_y * 1.4)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        hist_output_filename = f"{hist_output_folder}/{key}_comparison_hist.png"
        # Annotate yields for histogram plots (only for the plotted channel)
        ax = fig.axes[0]
        text_lines = []
        for cr_mode in sorted(cr_mode_data.keys()):
            vals = cr_mode_data[cr_mode]['unnorm_bins'].get(key, [])
            s = f"{float(np.sum(vals)):.1f}" if len(vals) > 0 else "0.0"
            text_lines.append(f"CR{cr_mode}: {s}")
        text_block = "\n".join(text_lines).rstrip()
        # place text left of the legend, right-aligned
        ax.text(0.78, 0.98, text_block, transform=ax.transAxes, fontsize=8, va='top', ha='right', family='monospace', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        ax.set_xlim(*resolve_x_limits(rivet_config_value))
        if rivet_config_value == "LHCb_default":
            ax.set_ylim(0, 0.4)
        fig.savefig(hist_output_filename, transparent=True)
        plt.close(fig)
        print(f"Plot saved: {hist_output_filename}")

    fig, ax = plt.subplots(figsize=(8, 6))
    plotted = False
    for cr_mode, cr_dict in cr_mode_data.items():
        normed = cr_dict["normalized"]
        total_data = normed.get('Total', None)
        if total_data:
            color = 'red' if cr_mode == '0' else 'blue'
            ax.stairs(
                np.array(total_data),
                np.array(edges),
                color=color,
                label=f"CRMode {cr_mode}",
                linewidth=2.2,
                fill=False,
            )
            plotted = True

    ax.set_xlabel(r'z(J/$\psi$)', fontsize=12)
    ax.set_ylabel('Normalized entries', fontsize=12)
    ax.set_title(f'{system_value}, {Energy_value}, ptHat {ptHat_value}, LDMEfac {ldme_value}, {event_value}, {rivet_config_value}\nAll Channels Comparison (Histogram)', fontsize=13)
    if plotted:
        ax.legend(loc='upper right', fontsize=10)
    # Dynamic y-axis for Total histogram: 1.4 * max(total) across CR modes
    max_y = 0.0
    for cr_mode, cr_dict in cr_mode_data.items():
        normed_local = cr_dict.get('normalized', {})
        total_data = normed_local.get('Total', None)
        if total_data:
            vals = np.array(total_data)
            if vals.size:
                max_y = max(max_y, float(np.max(vals)))
    if max_y <= 0.0:
        max_y = 1.0
    ax.set_ylim(0, max_y * 1.4)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    ax.set_xlim(*resolve_x_limits(rivet_config_value))
    if rivet_config_value == "LHCb_default":
        ax.set_ylim(0, 0.4)
    total_hist_output_filename = f"{hist_output_folder}/comparison_CRModes_hist.png"
    fig.savefig(total_hist_output_filename, transparent=True)
    plt.close(fig)
    print("\nComparison plot saved to '{}'".format(total_hist_output_filename))

# ====================================================================
#                           MAIN EXECUTION
# ====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compare CR modes across the J/psi z distribution using existing YODA files."
    )
    parser.add_argument("--energy", default="13TeV", help="Energy tag used in the input filename")
    parser.add_argument("--system", default="pp", help="System tag used in the input filename")
    parser.add_argument("--ptHat", default="25", help="ptHat tag used in the input filename")
    parser.add_argument("--ldme", default="1", help="LDMEfac tag used in the input filename")
    parser.add_argument("--event", default="10000000", help="Event tag used in the input filename")
    parser.add_argument("--rivet-config", default="LHCb_default", help="Rivet configuration tag used to pick the raw histogram path")
    parser.add_argument("--cr-modes", nargs="*", default=["0", "1"], help="CR modes to compare")
    parser.add_argument(
        "--input-template",
        default="/home/bayuadinp/Analyses/OniaShower/Production/rivetFile/StandAlone/{Energy_value}/{System_value}_OniaShower{channel}_{Energy_value}_ptHat{ptHat_value}_LDMEfac{LDME_value}_CR{CRMode}_{Event_value}.yoda",
        help="Filename template for the input YODA files",
    )
    parser.add_argument(
        "--output-folder",
        default=None,
        help="Optional output folder. Defaults to fig/CRMode_Comparisons/<system>_<energy>_LDMEfac<ldme>_ptHat<ptHat>_<event>",
    )
    parser.add_argument(
        "--normalize-output-folder",
        default=None,
        help="Optional folder for normalized YODA outputs. Defaults to Plotter/normalize.",
    )
    args = parser.parse_args()

    # 1. Define configuration variables here (consistent casing!)
    ENERGY = args.energy
    SYSTEM = args.system
    PTHAT = args.ptHat
    LDME = args.ldme
    EVENT = args.event
    RIVET_CONFIG = args.rivet_config

    CR_MODES = args.cr_modes
    # Update template to use the new variable names
    INPUT_FILES_TEMPLATE = args.input_template
    # INPUT_FILES_TEMPLATE = "StandAlone/{Energy_value}/{System_value}_OniaShower{channel}_{Energy_value}_ptHat{ptHat_value}_LDMEfac{LDME_value}_CR{CRMode}_{Event_value}_AllAnalyses.yoda"
    HIST_PATH = resolve_hist_path(RIVET_CONFIG)
    

    cr_mode_data = {}
    edges = []

    for cr_mode in CR_MODES:
        component_data = {}
        component_unnorm = {}
        found_any_channel = False
        
        # We check each channel independently
        for channel in ["Pwave", "Swave", "Singlet", "Octet"]:
            filename = INPUT_FILES_TEMPLATE.format(
                channel=channel, Energy_value=ENERGY, System_value=SYSTEM,
                ptHat_value=PTHAT, LDME_value=LDME, CRMode=cr_mode,
                Event_value=EVENT, RivetConfig_value=RIVET_CONFIG
            )

            if not os.path.exists(filename):
                print(f"Skipping: {channel} for CRMode {cr_mode} (File not found).")
                continue

            content = read_file_content(filename)
            values, scale = parse_yoda_values_and_scale(content, HIST_PATH)
            
            if values:
                # Capture edges from the first valid file found
                if not edges:
                    edges, _, _ = get_yoda_template_info(content, HIST_PATH)
                
                # Ensure data matches the expected bin count
                if edges and len(values) > (len(edges) - 1):
                    values = values[:len(edges) - 1]
                
                component_data[channel] = values
                found_any_channel = True

        if not found_any_channel:
            print(f"No data found for CRMode {cr_mode} at all.")
            continue

        # Fill missing channels with zeros so the calculation function doesn't crash
        # This allows "Total" to still be calculated from available parts
        num_bins = len(edges) - 1
        for channel in ["Pwave", "Swave", "Singlet", "Octet"]:
            if channel not in component_data:
                component_data[channel] = [0.0] * num_bins

        # Now calculate normalization normally
        norm_dict, norm_factors, unnorm_bins = calculate_normalized_data(
            component_data["Pwave"], component_data["Swave"],
            component_data["Singlet"], component_data["Octet"]
        )

        cr_mode_data[cr_mode] = {
            "normalized": norm_dict,
            "norm_factors": norm_factors,
            "unnorm_bins": unnorm_bins
        }

    normalize_out_folder = args.normalize_output_folder or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "normalize"
    )
    write_normalized_yoda_outputs(
        edges,
        cr_mode_data,
        normalize_out_folder,
        SYSTEM,
        ENERGY,
        LDME,
        PTHAT,
        EVENT,
        RIVET_CONFIG,
    )

    # Plotting logic remains the same
    out_folder = args.output_folder or f"fig/CRMode_Comparisons/{SYSTEM}_{ENERGY}_LDMEfac{LDME}_ptHat{PTHAT}_{EVENT}_rivet_config_{RIVET_CONFIG}"
    plot_comparison_for_cr_modes(edges, cr_mode_data, out_folder, SYSTEM, ENERGY, LDME, PTHAT, EVENT, RIVET_CONFIG)

if __name__ == "__main__":
    main()