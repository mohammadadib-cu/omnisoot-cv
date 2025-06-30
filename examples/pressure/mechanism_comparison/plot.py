from collections import defaultdict
from pathlib import Path
import os
import shutil

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
    "volume_fraction" : ["volume_fraction"], 
}

MW_carbon = 0.012; # kg/mol

def fetch_species_array(df, key):
    species_names = HEADER_DICT[key];
    for sp in species_names:
        if sp in df.columns:
            return df[sp].to_numpy();
    return None;

def save_plot(name, postprocess_dir, pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    # plt.tight_layout(pad = pad);
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

def use_tex():
    """Check if LaTeX is installed on the system."""
    if shutil.which('latex'):
        return True;
    else:
        return False;

def plot_nosoot():    
    mech_names = ["FFCM2", "Caltech", "KAUST", "ABF1bar", "ITV", "CRECK"];
    # Plot configs
    result_dir = f"results/";
    exp_dir = "data";
    postprocess_dir = f"post_process";
    
    linestyles = ["solid", "dashed", "dashdot", "dotted"];
    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][0];

    color_offset = 1
    step = 1;
    n_models = len(mech_names);
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)

    colors = [cmap(color_offset + i*step) for i in range(n_models)];
    color_black = (np.float64(0), np.float64(0), np.float64(0), np.float64(1))
    colors = [color_black] + colors[1:];
    
    exp_color = "#ff8833";
    area_alpha = 0.2;
    linewidth = 2;
    figsize = (6, 5.5);
    fontsize = 24;
    axis_fontsize = 26;
    plt.rcParams.update({"text.usetex": use_tex(), 'font.size': fontsize, 'figure.constrained_layout.use': True});
    
    labels = {
        "Caltech" : "Caltech",
        "KAUST" : "KAUST",
        "ABF1bar" : "ABF",
        "ITV" : "ITV",
        "CRECK" : "CRECK",
        "FFCM2" : "FFCM2",
        "CH4" : "$\\mathrm{CH_4}$", 
        "C2H4" : "$\\mathrm{C_2H_4}$",
        "C2H2" : "$\\mathrm{C_2H_2}$",
    }

    # Reading the simulation results
    results_df = defaultdict(
        dict
    );
    
    for mech_name in mech_names:
        file_path =  f"results/{mech_name}/sim_results.csv"
        df = pd.read_csv(file_path);
        results_df[mech_name] = df;
    

    # Reading the experimental data
    exp_df = pd.read_csv("data/species.csv");

    # ---------------------------------------------------------------------------------------------------------------
    # Ploting mole fraction of CH4
    # ---------------------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize);
    
    header = "CH4";
    scale_factor = 100;
    time_factor = 1000;
    
    # Experimental
    t_data_exp = exp_df["t [s]"].to_numpy();
    f_data_exp = exp_df["X_CH4 [-]"].to_numpy();
    f_sigma_exp = exp_df["sigmaX_CH4 [-]"].to_numpy();

    ax.plot(t_data_exp*time_factor, f_data_exp*scale_factor,
        linewidth = linewidth, linestyle = "solid", 
            label = "Data", color = exp_color);

    ax.fill_between(t_data_exp*time_factor, (f_data_exp - f_sigma_exp)*scale_factor, (f_data_exp + f_sigma_exp)*scale_factor,
                color = exp_color, alpha = area_alpha)
    
    # Numerical
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data = fetch_species_array(df, header);
        linestyle = "dashed" if mech_name == "FFCM2" else "solid"
        ax.plot(t_data*time_factor, f_data*scale_factor,
            linewidth = linewidth, linestyle = linestyle,
                label = labels[mech_name], color = colors[mech_index]);
    
    #Y axis
    ax.set_ylabel(f"Mole fraction of " + labels[header] + scale_to_label(scale_factor), {'fontsize' : axis_fontsize});
    ax.set_ylim([0, 5.1]);
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.2))

    # X axis
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    ax.set_xlim([-0.03, 2.03]);
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.legend(loc = (0.55, 0.35), ncol = 1, fontsize = 19, columnspacing = 0.1,
                    framealpha=0);

    ax.text(0.15, 0.78, "(a)", transform=ax.transAxes, fontsize = 32)

    save_plot(header, postprocess_dir = postprocess_dir, pad = 0.4)


    # ---------------------------------------------------------------------------------------------------------------
    # Ploting the carbon mass fraction of C2H4
    # ---------------------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize);
    header = "C2H4";
    scale_factor = 1000;
    time_factor = 1000;

    mech_index = 0;
    mech_name = mech_names[mech_index]

    # Experimental
    t_data_exp = exp_df["t [s]"].to_numpy();
    f_data_exp = exp_df["X_C2H4 [-]"].to_numpy();
    f_sigma_exp = exp_df["sigmaX_C2H4 [-]"].to_numpy();

    ax.plot(t_data_exp*time_factor, f_data_exp*scale_factor,
        linewidth = linewidth, linestyle = "solid", 
            label = "Data", color = exp_color);

    ax.fill_between(t_data_exp*time_factor, (f_data_exp - f_sigma_exp)*scale_factor, (f_data_exp + f_sigma_exp)*scale_factor,
                color = exp_color, alpha = area_alpha)

    # Numerical
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data =fetch_species_array(df, header);
        linestyle = "dashed" if mech_name == "FFCM2" else "solid"
        ax.plot(t_data*time_factor, f_data*scale_factor,
            linewidth = linewidth, linestyle = linestyle,
                label = labels[mech_name], color = colors[mech_index]);

    #Y axis
    ax.set_ylabel(f"Mole fraction of " + labels[header] + scale_to_label(scale_factor), {'fontsize' : axis_fontsize});
    ax.set_ylim([0, 8.5]);
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.4))

    # X axis
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    ax.set_xlim([-0.03, 2.03]);
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.text(0.1, 0.85, "(b)", transform=ax.transAxes, fontsize = 32);

    #### Zoomed in ------------------------------

    left, bottom, width, height = [0.4, 0.42, 0.43, 0.4]
    ax2 = fig.add_axes([left, bottom, width, height])

    # Experimental
    t_data_exp = exp_df["t [s]"].to_numpy();
    f_data_exp = exp_df["X_C2H4 [-]"].to_numpy();
    f_sigma_exp = exp_df["sigmaX_C2H4 [-]"].to_numpy();

    ax2.plot(t_data_exp*time_factor, f_data_exp*scale_factor,
        linewidth = linewidth, linestyle = "solid", 
            label = "Data", color = exp_color);

    ax2.fill_between(t_data_exp*time_factor, (f_data_exp - f_sigma_exp)*scale_factor, (f_data_exp + f_sigma_exp)*scale_factor,
                color = exp_color, alpha = area_alpha)

    # Numerical
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data =fetch_species_array(df, header);
        linestyle = "dashed" if mech_name == "FFCM2" else "solid"
        ax2.plot(t_data*time_factor, f_data*scale_factor,
            linewidth = linewidth, linestyle = linestyle,
                label = labels[mech_name], color = colors[mech_index]);

    #Y axis
    ax2.set_ylim([0, 8.5]);
    ax2.yaxis.set_major_locator(MultipleLocator(2))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.4))
    ax2.set_yticklabels(ax2.get_yticklabels(), fontsize = 14)

    # X axis
    ax2.set_xlim([-0.005, 0.2]);
    ax2.xaxis.set_major_locator(MultipleLocator(0.1))
    ax2.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax2.set_xticklabels(ax2.get_xticklabels(), fontsize = 14)

    save_plot(header, postprocess_dir = postprocess_dir, pad = 0.4)

    # ---------------------------------------------------------------------------------------------------------------
    # Ploting the carbon mass fraction of C2H2
    # ---------------------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize);
    header = "C2H2";
    scale_factor = 100;
    time_factor = 1000;
    # Experimental
    t_data_exp = exp_df["t [s]"].to_numpy();
    f_data_exp = exp_df["X_C2H2 [-]"].to_numpy();
    f_sigma_exp = exp_df["sigmaX_C2H2 [-]"].to_numpy();

    ax.plot(t_data_exp*time_factor, f_data_exp*scale_factor,
        linewidth = linewidth, linestyle = "solid", 
            label = "Data", color = exp_color);

    ax.fill_between(t_data_exp*time_factor, (f_data_exp - f_sigma_exp)*scale_factor, (f_data_exp + f_sigma_exp)*scale_factor,
                color = exp_color, alpha = area_alpha)

    # Numerical
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data =fetch_species_array(df, header);
        linestyle = "dashed" if mech_name == "FFCM2" else "solid"
        ax.plot(t_data*time_factor, f_data*scale_factor,
            linewidth = linewidth, linestyle = linestyle,
                label = labels[mech_name], color = colors[mech_index]);

    #Y axis
    ax.set_ylabel(f"Mole fraction of " + labels[header] + scale_to_label(scale_factor), {'fontsize' : axis_fontsize});
    ax.set_ylim([0, 2.15]);
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    # X axis
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    ax.set_xlim([-0.03, 2.03]);
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.text(0.06, 0.87, "(c)", transform=ax.transAxes, fontsize = 32);
    
    save_plot(header, postprocess_dir = postprocess_dir, pad = 0.4)


def plot_maxsoot():
    # Plot configs
    exp_color = "#ff8833";
    area_alpha = 0.2
    linewidth = 2;
    figsize = (6, 6.5);
    fontsize = 24;
    axis_fontsize = 26;
    plt.rcParams.update({"text.usetex": use_tex(), 'font.size': fontsize, 'figure.constrained_layout.use': True});
    
    result_dir = f"results_Max/";
    exp_dir = "data";
    postprocess_dir = f"post_process";
    labels = {
        "Caltech" : "Caltech",
        "KAUST" : "KAUST",
        "ABF1bar" : "ABF",
        "ITV" : "ITV",
        "CRECK" : "CRECK",
    };
    mech_names = ["Caltech", "KAUST", "ABF1bar", "ITV", "CRECK"];
    # Reading the simulation results
    results_df = defaultdict(
        dict
    );

    # Reading the experimental data
    exp_df = pd.read_csv(f"{exp_dir}/soot.csv");

    for mech_name in mech_names:
        file_path =  f"{result_dir}/{mech_name}/IrreversibleDimerization//T2203K.csv"
        df = pd.read_csv(file_path);
        results_df[mech_name] = df;


    fig, ax = plt.subplots(figsize=figsize);

    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][0];
    color_offset = 1
    step = 1;
    n_models = len(mech_names)+1;
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)

    colors = [cmap(color_offset + i*step) for i in range(n_models)];
    colors = colors[1:];

    header = "volume_fraction";

    Em_s = [0.174, 0.37]

    scale_factor = 1e6;
    time_factor = 1000;

    mech_index = 0;
    mech_name = mech_names[mech_index];

    # Experimental
    t_data_exp = exp_df["t [s]"].to_numpy();
    f_rel_data_exp = exp_df["fv_rel_633nm"].to_numpy();
    f_max_data_exp, f_min_data_exp = f_rel_data_exp / Em_s[0], f_rel_data_exp / Em_s[1];
    f_mean_data_exp = (f_max_data_exp + f_min_data_exp) / 2.0;

    ax.fill_between(t_data_exp*time_factor, f_min_data_exp*scale_factor, f_max_data_exp*scale_factor,
                color = exp_color, alpha = area_alpha)
    ax.plot(t_data_exp*time_factor, f_mean_data_exp*scale_factor,
        linewidth = linewidth, linestyle = "solid", 
            label = "Data", color = exp_color);
        

    # Numerical
    for mech_index, mech_name in enumerate(mech_names):
        df = results_df[mech_name];
        t_data = df["t[s]"];
        f_data = fetch_species_array(df, header);
        if not (f_data is None):
            ax.plot(t_data*time_factor, f_data*scale_factor,
                linewidth = linewidth, linestyle = "solid",
                    label = labels[mech_name], color = colors[mech_index]);




    #Y axis
    ax.set_ylabel("Soot volume fraction, $f_v$ [ppm]", {'fontsize' : axis_fontsize});
    ax.set_ylim([0, 2]);
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    # X axis
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    ax.set_xlim([-0.05, 2.01]);
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax.legend(loc = (0.1, 1), ncol = 2, fontsize = 19, columnspacing = 1.1,
                    framealpha=0);

    save_plot("soot_fv_max", postprocess_dir = postprocess_dir, pad = 0.4)


if __name__ == "__main__":
    plot_nosoot();
    plot_maxsoot();

