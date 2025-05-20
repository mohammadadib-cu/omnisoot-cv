import os
from collections import defaultdict
from pathlib import Path
import string
import shutil

import cantera as ct
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, LogLocator)
import numpy as np
import pandas as pd

def save_plot(name, postprocess_dir, save_pdf = True, white_back = "None", pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    plt.tight_layout(pad = pad);
    facecolor = white_back;
    plt.savefig(f"{postprocess_dir}/{name}.svg", facecolor=facecolor);
    if save_pdf:
        filename = f"{postprocess_dir}/{name}.pdf";
        plt.savefig(filename);
        
def plot():
    # Check if Latex is installed
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

    font_family = 'Times Modern Roman';
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
    n_secs = 60;

    # Reading results
    results_df = {};
    for particle_dynamics_model in particle_dynamics_models:
        results_df[particle_dynamics_model] = pd.read_csv(f"{result_dir}/{particle_dynamics_model}.csv");


    # Reading data
    data_dir = "data";
    df_DEM = pd.read_csv(f"{data_dir}/Armen.csv")
    df_sigmag = pd.read_csv(f"{data_dir}/sigmag_dm.csv")
    df_eta_psi = pd.read_csv(f"{data_dir}/DEM_eta_psi.csv")

    # ----------------------------------------------------------------------------------------------------------------
    # Plotting N_{agg} & N_{pri}
    # ----------------------------------------------------------------------------------------------------------------
    plt.rcParams.update({"text.usetex": use_tex, 'font.family':'Latin Modern Roman', 'font.size': fontsize});
    fig, ax = plt.subplots(figsize=figsize);
    
    for particle_dynamics_model_index, particle_dynamics_model in enumerate(particle_dynamics_models):
        t_data = results_df[particle_dynamics_model]["t[s]"] * 1000;
        z_data = results_df[particle_dynamics_model]["N_pri[#/m3]"];
        ax.plot(t_data, z_data,
            color=colors[particle_dynamics_model_index], linewidth = linewidth+particle_dynamics_model_index*0.5, label = " ",
                linestyle = linestyles[1]);
    
    # Placeholder
    ax.plot([0], [0], label = " ", alpha = 0);
    
    
    for particle_dynamics_model_index, particle_dynamics_model in enumerate(particle_dynamics_models):
        t_data = results_df[particle_dynamics_model]["t[s]"] * 1000;
        z_data = results_df[particle_dynamics_model]["N_agg[#/m3]"];
        ax.plot(t_data, z_data,
            color=colors[particle_dynamics_model_index], linewidth = linewidth+particle_dynamics_model_index*0.5, label = labels[particle_dynamics_model],
                linestyle = linestyles[0]);
    
    t_data = df_DEM['time'][0:-1:200]*1000;
    z_data = df_DEM['Np'][0:-1:200];
    ax.plot(t_data, z_data,
        color="#1d6e1f", markeredgewidth = markeredgewidth, label = labels["DEM"], linestyle = "none", 
            marker = 'o', markersize = markersize, markerfacecolor="#ffffff");
    
    ax.set_yscale("log");
    ax.set_ylim([1e15, 6e18]);
    ax.set_ylabel("$N\\;\\mathrm{[1/m^3]}$", {'fontsize' : axis_fontsize});
    
    ax.set_xlim([-0.6, 35.01]);
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    
    xticks = list(range(0, 36, 5));
    ax.set_xticks(xticks);
    ax.xaxis.set_minor_locator(MultipleLocator(1));
    
    ax.text(0.47, 0.67, "${N_{agg}}$", transform=ax.transAxes, fontsize = legend_fontsize+4)
    ax.text(0.28, 0.67, "${N_{pri}}$", transform=ax.transAxes, fontsize = legend_fontsize+4)
    legend = ax.legend(loc = (0.25, 0.35), ncol=2, fontsize = legend_fontsize+3, framealpha=0, columnspacing=0.3);
    
    ax.text(0.85, 0.75, "(a)", transform=ax.transAxes, fontsize = 28)
    
    save_plot(f"N_agg_pri", postprocess_dir = postprocess_dir)
    # ----------------------------------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------------------------------
    # Plotting d_m & d_g
    # ----------------------------------------------------------------------------------------------------------------    
    plt.rcParams.update({"text.usetex": use_tex, 'font.family': font_family, 'font.size': fontsize});
    fig, ax = plt.subplots(figsize=figsize);
    
    # Mobility
    for particle_dynamics_model_index, particle_dynamics_model in enumerate(particle_dynamics_models):
        df = results_df[particle_dynamics_model]
        t_data = df["t[s]"] * 1000;
        z_data = df["d_m[nm]"];
        ax.plot(t_data, z_data,
            color=colors[particle_dynamics_model_index], linewidth = linewidth+0.5*particle_dynamics_model_index, label = " ",
                linestyle = linestyles[0]);
    
    t_data = df_DEM['time'][0:-1:300]*1000;
    z_data = df_DEM['dm '][0:-1:300]*1e9;
    ax.plot(t_data, z_data,
        color="#308A33", linewidth = linewidth, label = " ", linestyle = "none", 
            marker = 'o', markersize = markersize, markerfacecolor="#ffffff",
           markeredgewidth = markeredgewidth,
           );
    
    # Gyration
    for particle_dynamics_model_index, particle_dynamics_model in enumerate(particle_dynamics_models):
        df = results_df[particle_dynamics_model]
        t_data = df["t[s]"]*1000;
        z_data = df["d_g[nm]"];
        ax.plot(t_data, z_data,
            color=colors[particle_dynamics_model_index], linewidth = linewidth+0.5*particle_dynamics_model_index, 
                label = labels[particle_dynamics_model], linestyle = linestyles[1],
               );
    
    t_data = df_DEM['time'][0:-1:300]*1000;
    z_data = df_DEM['dg'][0:-1:300]*1e9;
    ax.plot(t_data, z_data,
        color="#a3009e", linewidth = linewidth, label = labels["DEM"], linestyle = "none", 
            marker = 's', markersize = markersize,
           markeredgewidth = 1, markerfacecolor="#ffffff");
    
    ax.set_ylim([0, 75]);
    ax.set_ylabel("${d}$ [nm]", {'fontsize' : axis_fontsize});
    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(4));
    
    ax.set_xlim([-0.15, 35.01]);
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    xticks = list(range(0, 36, 5));
    ax.set_xticks(xticks);
    ax.xaxis.set_minor_locator(MultipleLocator(1));
    
    x_offset = 0.05;
    y_offset = 0.2
    ax.text(0.52 + x_offset, 0.33, "${d_{g}}$", transform=ax.transAxes, fontsize = legend_fontsize + 5)
    ax.text(0.39 + x_offset, 0.33, "${d_{m}}$", transform=ax.transAxes, fontsize = legend_fontsize + 5)
    legend = ax.legend(loc = (0.35 + x_offset, 0.02), ncol=2, columnspacing = 0.1, 
                       fontsize = legend_fontsize+3, framealpha=0, handletextpad = 0.6, handlelength = 1.6);
    
    ax.text(0.05, 0.85, "(b)", transform=ax.transAxes, fontsize = 28)
    save_plot(f"d_mg", postprocess_dir = postprocess_dir)
    # ----------------------------------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------------------------------
    # Plotting sigma_{mg}
    # ----------------------------------------------------------------------------------------------------------------   
    plt.rcParams.update({"text.usetex": use_tex,'font.family': font_family, 'font.size': fontsize});
    fig, ax = plt.subplots(figsize = figsize);
    
    particle_dynamics_model = "Sectional"
    df = results_df[particle_dynamics_model]
    
    t_data = df["t[s]"]*1000;
    z_data = df["sigma_mg"];
    ax.plot(t_data, z_data,
        color= "#e30000", linewidth = linewidth, label = labels[particle_dynamics_model], linestyle = linestyles[0]);
    
    
    t_data = df_sigmag['t']*1000;
    z_data = df_sigmag['sigmag_dm'];
    ax.plot(t_data, z_data,
        color="#308A33", linewidth = linewidth, label = labels["DEM"], linestyle = "none", 
            marker = 'o', markersize = markersize-1, markerfacecolor="none",
           markeredgewidth = 1.5,
           alpha = 1);
    
    t_data = [8.5, 13];
    z_data = [2.03, 2.03];
    ax.plot(t_data, z_data,
        color="#b50082", linewidth = linewidth, linestyle = "dashed")
    
    ax.text(8.3, 2.2, "${\\sigma_{g, fm}}$=2.03", fontsize = 22, color="#b50082")
    
    ax.set_ylim([0.5, 2.5]);
    ax.yaxis.set_major_locator(MultipleLocator(0.5));
    ax.yaxis.set_minor_locator(MultipleLocator(0.1));
    ax.set_ylabel("Mobility diameter\nstandard deviation, ${\\sigma_{g}}$", {'fontsize' : axis_fontsize});
    

    ax.set_xlim([-0.2, 12.1]);
    ax.set_xlabel("t [ms]", {'fontsize' : axis_fontsize});
    xticks = list(range(0, 13, 2));
    ax.set_xticks(xticks);
    ax.xaxis.set_minor_locator(MultipleLocator(0.4));
    
    legend = ax.legend(loc = (0.15, 0.1), fontsize = legend_fontsize, framealpha=0);
    
    ax.text(0.05, 0.85, "(a)", transform=ax.transAxes, fontsize = 28)
    
    save_plot("sigmag", postprocess_dir = postprocess_dir)
    # ----------------------------------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------------------------------
    # Plotting sigma_{mg}
    # ---------------------------------------------------------------------------------------------------------------- 
    plt.rcParams.update({"text.usetex": use_tex, 'font.family': font_family, 'font.size': fontsize});
    fig, ax = plt.subplots(figsize = figsize);
    
    particle_dynamics_model = "Sectional"
    df = results_df[particle_dynamics_model];
    
    # Sampling times for SPSD
    t_samples = [1e-3, 4e-3, 22e-3, 450e-3, 680e-3];
    
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
    
    t_data = df_eta_psi["eta"]
    z_data = df_eta_psi["psi"]
    ax.plot(t_data, z_data, color="#000000", label = "DEM", linestyle = "none", 
            marker = 'o', markersize = markersize-1, markerfacecolor="none",
           markeredgewidth = 1);
    
    ax.set_yscale("log");
    ax.set_ylim([1e-6, 1]);
    ax.set_ylabel("Normalized concentration, $\\mathrm{\\Psi \\eta}$", {'fontsize' : axis_fontsize});
    
    
    ax.set_xscale("log");
    ax.set_xlabel("Dimensionless volume, $\\mathrm{\\eta}$", {'fontsize' : axis_fontsize});
    ax.set_xlim([1e-3, 100]);
        
    legend = ax.legend(loc = (0.3, 0.05), fontsize = legend_fontsize+3, framealpha=0);
    ax.text(0.05, 0.85, "(b)", transform=ax.transAxes, fontsize = 28);
    
    save_plot(f"PSD", postprocess_dir = postprocess_dir)