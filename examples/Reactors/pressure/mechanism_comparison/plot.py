from collections import defaultdict
from pathlib import Path
import os

import cantera as ct
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import pandas as pd

MECH_DICT = {
    "Caltech":"./data/Caltech.yaml",
    "KAUST": "./data/KAUST.yaml",
    "ABF1bar" : "./data/abf_1bar.yaml",
    "CRECK" : "./data/CRECK_no_BINS.yaml",
    "ITV" : "./data/ITV_PAH.yaml",
    "NUIG" : "./data/NUIG.yaml",
    "FFCM2" : "./data/FFCM2.yaml",
}

HEADER_DICT = {
    "CH4" : ["CH4"],
    "C2H4" : ["C2H4"],
    "C2H2" : ["C2H2"],
    "A1" : ["A1", "C6H6"],
    "A2" : ["A2", "C10H8"],
    "A2" : ["A2", "C10H8"],
    "A3" : ["A3", "C14H10", "A3XC14H10"],
    "A2R5" : ["A2R5", "C12H8"], 
    "A4" : ["A4", "C16H10", "A4XC16H10"], 
    "A3R5" : ["A3R5"], 
    "A4R5" : ["A4R5"], 
}

MW_carbon = 0.012; # kg/mol

def fetch_species_array(df, key):
    species_names = HEADER_DICT[key];
    for sp in species_names:
        if sp in df.columns:
            return df[sp].to_numpy();
    return None;

def save_plot(name, postprocess_dir, save_pdf = True, white_back = "None", in_paper = False, pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    # plt.tight_layout(pad = pad);
    facecolor = white_back;
    plt.savefig(f"{postprocess_dir}/{name}.svg", facecolor=facecolor);
    if save_pdf:
        filename = f"{postprocess_dir}/{name}.pdf";
        plt.savefig(filename);


def scale_to_label(scale_factor):
    exponent = int(np.log10(scale_factor));
    if exponent == 0:
        return "";
    elif exponent == 1:
        return "$\\times$10"
    return (
        f"$\\times10^{exponent}$"
    )


def plot():
    mech_names = ["NUIG", "FFCM2", "Caltech", "KAUST", "ABF1bar", "ITV", "CRECK"];

    # Plot configs
    result_dir = f"results/";
    exp_dir = "data";
    postprocess_dir = f"post_process";
    
    linestyles = ["solid", "dashed", "dashdot", "dotted"];
    
    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][0];
    
    color_offset = 1
    step = 1;
    n_models = len(mech_names) - 1;
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)
    colors = [cmap(color_offset + i*step) for i in range(n_models)];
    color_black = (np.float64(0), np.float64(0), np.float64(0), np.float64(1))
    colors = colors[:1] + [color_black] + colors[1:];
    
    linewidth = 2;
    markersize = 10;
    ticklength = 5;
    tickwidth = 0.75;
    figsize = (6, 5.5);
    fontsize = 24;
    legend_fontsize = 18;
    axis_fontsize = 26;
    letter_fontsize = 32;
    
    labels = {
        "Caltech" : "Caltech",
        "KAUST" : "KAUST",
        "ABF1bar" : "ABF",
        "ITV" : "ITV",
        "CRECK" : "CRECK",
        "NUIG" : "NUIG",
        "FFCM2" : "FFCM2",
        "CH4" : "$\\mathrm{CH_4}$", 
        "C2H4" : "$\\mathrm{C_2H_4}$",
        "C2H2" : "$\\mathrm{C_2H_2}$",
        "A1" : "A1",
        "A2" : "A2",
        "A3" : "A3",
        "A4" : "A4",
        "A2R5" : "A2R5",
        "A3R5" : "A3R5",
        "A4R5" : "A4R5",
    }

    # Reading the simulation results
    results_df = defaultdict(
        dict
    );
    
    for mech_name in mech_names:
        file_path =  f"results/{mech_name}/sim_results.csv"
        df = pd.read_csv(file_path);
        results_df[mech_name] = df;
    
    # Calculating the initial carbon mass fraction
    df = results_df["Caltech"];
    sp_index = list(df.columns).index("N2");
    x_columns = list(df.columns[sp_index:]);
    gas = ct.Solution(MECH_DICT["Caltech"]);
    lookup_index = 0;
    gas.TPX = df.iloc[lookup_index]["T[K]"].item(), df.iloc[lookup_index]["P[Pa]"].item(), df.iloc[lookup_index][x_columns].to_numpy();
    mfc_c_0 = 0.0;
    n_c_dict = {sp:gas.n_atoms(sp, "C") for sp in gas.species_names}
    for sp_index, sp in enumerate(gas.species_names):
        mfc_c_0 += gas.n_atoms(sp, "C") * MW_carbon * gas.X[sp_index].item() / (gas.mean_molecular_weight / 1000)

    # ---------------------------------------------------------------------------------------------------------------
    # Ploting the carbon mass fraction of CH4
    # ---------------------------------------------------------------------------------------------------------------
    plt.rcParams.update({"text.usetex": True, 'font.family':'Latin Modern Roman', 'font.size': fontsize, 'figure.constrained_layout.use': True});
    fig, ax = plt.subplots(figsize=figsize);
    
    header = "CH4";
    scale_factor = 1;
    time_factor = 1000;
    
    mech_index = 0;
    mech_name = mech_names[mech_index]
    
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data = n_c_dict[header] * fetch_species_array(df, header) * MW_carbon / df["mean_MW[kg/mol]"] / mfc_c_0;
        ax.plot(t_data*time_factor, f_data*scale_factor,
               linewidth = linewidth, linestyle = linestyles[0],
                label = labels[mech_name], color = colors[mech_index]);
    
    #Y axis
    ax.set_ylabel(f"Carbon mass fraction of " + labels[header] + scale_to_label(scale_factor), {'fontsize' : axis_fontsize});
    ax.set_ylim([0, 0.8]);
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.04))
    
    # X axis
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    ax.set_xlim([-0.1, 5.1]);
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    legend = ax.legend(loc = (0.61, 0.4), ncol = 1, fontsize = 18, columnspacing = 0.1,
                    framealpha=0);
    
    left, bottom, width, height = [0.28, 0.52, 0.35, 0.4]
    ax2 = fig.add_axes([left, bottom, width, height])
    
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data = n_c_dict[header] * fetch_species_array(df, header) * MW_carbon / df["mean_MW[kg/mol]"] / mfc_c_0;
        if not (f_data is None):
            ax2.plot(t_data*time_factor, f_data*scale_factor,
                   linewidth = linewidth, linestyle = linestyles[0],
                    label = labels[mech_name], color = colors[mech_index]);
    
    ax2.tick_params(axis='both', labelsize = 14)
    #Y axis
    ax2.set_ylim([0.2, 1]);
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.04));
    
    # X axis
    ax2.set_xlim([1e-3, 1e0]);
    ax2.set_xscale("log");
    #ax2.set_xticklabels(ax2.get_xticklabels(), fontsize = 14);
    
    ax2.text(0.81, 0.28, "(a)", transform=ax.transAxes, fontsize = 32)
    
    save_plot(header, postprocess_dir = postprocess_dir, pad = 0.4)


    # ---------------------------------------------------------------------------------------------------------------
    # Ploting the carbon mass fraction of C2H2
    # ---------------------------------------------------------------------------------------------------------------
    plt.rcParams.update({"text.usetex": True, 'font.family':'Latin Modern Roman', 'font.size': fontsize, 'figure.constrained_layout.use': True});
    fig, ax = plt.subplots(figsize=figsize);

    header = "C2H2";
    
    scale_factor = 1;
    time_factor = 1000;
    
    mech_index = 0;
    mech_name = mech_names[mech_index]
    
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data = n_c_dict[header] * fetch_species_array(df, header) * MW_carbon / df["mean_MW[kg/mol]"] / mfc_c_0;
        ax.plot(t_data*time_factor, f_data*scale_factor,
               linewidth = linewidth, linestyle = linestyles[0],
                label = labels[mech_name], color = colors[mech_index]);
    
    #Y axis
    ax.set_ylabel(f"Carbon mass fraction of " + labels[header] + scale_to_label(scale_factor), {'fontsize' : axis_fontsize});
    ax.set_ylim([0, 0.85]);
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.04))
    
    # X axis
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    ax.set_xlim([-0.1, 5.1]);
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    
    ax.text(0.8, 0.45, "(b)", transform=ax.transAxes, fontsize = 32);
    
    
    left, bottom, width, height = [0.34, 0.24, 0.35, 0.35]
    ax2 = fig.add_axes([left, bottom, width, height])
    
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data = n_c_dict[header] * fetch_species_array(df, header) * MW_carbon / df["mean_MW[kg/mol]"] / mfc_c_0;
        if not (f_data is None):
            ax2.plot(t_data*time_factor, f_data*scale_factor,
                   linewidth = linewidth, linestyle = linestyles[0],
                    label = labels[mech_name], color = colors[mech_index]);
    
    ax2.tick_params(axis='both', labelsize = 14)
    #Y axis
    ax2.set_ylim([0, 0.7]);
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.04))
    
    # X axis
    ax2.set_xlim([1e-3, 1e-0]);
    ax2.set_xscale("log");
    
    save_plot(header, postprocess_dir = postprocess_dir, pad = 0.4)

    # ---------------------------------------------------------------------------------------------------------------
    # Ploting the carbon mass fraction of A2R5
    # ---------------------------------------------------------------------------------------------------------------
    plt.rcParams.update({"text.usetex": True, 'font.family':'Latin Modern Roman', 'font.size': fontsize, 'figure.constrained_layout.use': True});
    fig, ax = plt.subplots(figsize=figsize);
    
    header = "A2R5";
    
    scale_factor = 1e2;
    time_factor = 1000;
    
    mech_index = 0;
    mech_name = mech_names[mech_index]
    
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data = fetch_species_array(df, header);
        if not (f_data is None):
            f_data = n_c_dict[header] * fetch_species_array(df, header) * MW_carbon / df["mean_MW[kg/mol]"] / mfc_c_0;
            ax.plot(t_data*time_factor, f_data*scale_factor,
                   linewidth = linewidth, linestyle = linestyles[0],
                    label = labels[mech_name], color = colors[mech_index]);
    
    
    #Y axis
    ax.set_ylabel(f"Carbon mass fraction of " + labels[header] + scale_to_label(scale_factor), {'fontsize' : axis_fontsize-3});
    ax.set_ylim([0, 5]);
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.2))
    
    # X axis
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    ax.set_xlim([-0.1, 5.1]);
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    
    
    legend = ax.legend(loc = (0.03, 0.72), ncol = 2, fontsize = 18, columnspacing = 0.5,
                    framealpha=0);
    
    ax.text(0.82, 0.84, "(c)", transform=ax.transAxes, fontsize = 32);
    
    save_plot(header, postprocess_dir = postprocess_dir, pad = 0.4)



if __name__ == "__main__":
    plot();

