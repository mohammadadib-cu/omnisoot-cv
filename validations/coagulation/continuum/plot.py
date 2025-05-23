import os
from collections import defaultdict
from pathlib import Path
import string
import shutil
import json

import cantera as ct
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, LogLocator)
import numpy as np
import pandas as pd

def save_plot(name, postprocess_dir = "", save_pdf = True, white_back = "None", pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    plt.tight_layout(pad = pad);
    facecolor = white_back;
    plt.savefig(f"{postprocess_dir}/{name}.svg", facecolor=facecolor);
    if save_pdf:
        filename = f"{postprocess_dir}/{name}.pdf";
        plt.savefig(filename);

def plot():
    if shutil.which('latex'):
        use_tex = True;
    else:
        use_tex = False;


    # Plot configs
    postprocess_dir = "post_process";
    result_dir = "results";

    linestyles = ["solid", "dashed", "dashdot", "dotted"];

    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][1];

    color_offset = 1
    step = 1;
    n_models = 2;
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)
    colors = [cmap(color_offset + i*step) for i in range(n_models)];

    linewidth = 1.5;
    markersize = 12;
    markeredgewidth = 2
    ticklength = 5;
    tickwidth = 0.75;
    figsize = (6.3, 5.2);
    fontsize = 22;
    legend_fontsize = 18;
    axis_fontsize = 24;
    letter_fontsize = 36;

    labels = {
        "Monodisperse" : "MPBM", 
        "Sectional" : "SPBM",
        "DEM" : "DEM"
    }

    particle_dynamics_models = ["Monodisperse", "Sectional"];
    n_secs = 80;
    sec_vars = ["eta", "psi"];
    sec_headers = {var: [f"{var}_{sec}" for sec in range(n_secs)] for var in sec_vars}

    # Reading the results
    results_df = {};
    for particle_dynamics_model in particle_dynamics_models:
        results_df[particle_dynamics_model] = pd.read_csv(f"{result_dir}/{particle_dynamics_model}.csv")

    # Reading the DEM results
    data_dir = "data";
    df_eta_psi = pd.read_csv(f"{data_dir}/DEM_eta_psi.csv");

    # Self-Preserving Size Distribution
    plt.rcParams.update({"text.usetex": use_tex,'font.family':'Latin Modern Roman', 'font.size': fontsize});
    fig, ax = plt.subplots(figsize=figsize);

    particle_dynamics_model = "Sectional"
    df = results_df[particle_dynamics_model];

    # Sampling times for SPSD
    t_samples = [20e-3, 80e-3, 200e-3, 500e-3, 900e-3];

    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][0];
    color_offset = 2
    step = 1;
    n_models = len(t_samples);
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)
    spsd_colors = [cmap(color_offset + i*step) for i in range(n_models)];

    for t_index, t_sample in enumerate(t_samples):
        eta_array = np.array([np.interp(t_sample, df["t[s]"], df[f"eta_{i}"]) for i in range(n_secs)])
        psi_array = np.array([np.interp(t_sample, df["t[s]"], df[f"psi_{i}"]) for i in range(n_secs)])
        label = f"t={(t_sample*1000):0.0f} ms";
        ax.plot(eta_array, eta_array * psi_array, label = label,
            color=spsd_colors[t_index], linewidth = linewidth, 
                linestyle = "solid");

    eta_exp = df_eta_psi["eta"]
    psi_exp = df_eta_psi["psi"]
    ax.plot(eta_exp, eta_exp * psi_exp, color="#000000", label = "DEM", linestyle = "none", 
            marker = 'o', markersize = markersize-1, markerfacecolor="none",
        markeredgewidth = 1);

    ax.set_yscale("log");
    ax.set_ylim([1e-5, 1]);
    ax.set_ylabel("Normalized concentration, $\\mathrm{\\Psi \\eta}$", {'fontsize' : axis_fontsize});


    ax.set_xscale("log");
    ax.set_xlabel("Dimensionless volume, $\\mathrm{\\eta}$", {'fontsize' : axis_fontsize});
    ax.set_xlim([1e-3, 100]);

    legend = ax.legend(loc = (0.33, 0.05), fontsize = legend_fontsize+3, framealpha=0, handlelength = 1.5);

    ax.text(0.05, 0.85, "(b)", transform=ax.transAxes, fontsize = 28);

    save_plot(f"PSD", postprocess_dir = postprocess_dir)