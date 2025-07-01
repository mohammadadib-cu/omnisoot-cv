from collections import defaultdict
from pathlib import Path
import os
import shutil
import string

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, LogLocator)
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

Av = 6.0221408e+23; # 1/mol

def save_plot(name, postprocess_dir, pad = 0.5):
    if not os.path.exists(postprocess_dir):
        os.makedirs(postprocess_dir);
    # plt.tight_layout(pad = pad);
    filename = f"{postprocess_dir}/{name}.pdf";
    plt.savefig(filename);

def use_tex():
    """Check if LaTeX is installed on the system."""
    if shutil.which('latex'):
        return True;
    else:
        return False;

def make_title(case):
    return (
        f"Q={case[:-6]} L/min"
    )


def plot():

    # Plot configurations
    postprocess_dir = "post_process";
    linestyles = ["solid", "dashed", "dashdot", "dotted"];
    cmap_name = ['hsv', "YlOrRd", "YlGnBu", "YlGn"][0];

    color_offset = 1
    step = 1;
    n_models = 4;
    cmap = plt.get_cmap(cmap_name, n_models*step + color_offset)

    colors = [cmap(color_offset + i*step) for i in range(n_models)];

    linewidth = 2;
    figsize = (6.1, 5.3);
    fontsize = 24;
    axis_fontsize = 26;
    letter_fontsize = 26;
    legend_fontsize = 20;
    title_fontsize = 26;

    plt.rcParams.update({"text.usetex": use_tex(), 'font.size': fontsize, 'figure.constrained_layout.use': True});

    labels = {
        "DimerCoalescence" : "Dim.Coal.",
        "ReactiveDimerization":"Reac.Dim.",
        "EBridgeModified" : "EBri.Mod.",
        "IrreversibleDimerization" : "Irrev.Dim.",
        "data" : "Data",
        "d_p" : "$d_p$",
        "d_m" : "$d_m$",
    };

    gas_density = 1.13153; # kg/m3
    cases = ["8.0lmin", "11.0lmin", "12.0lmin"];
    PAH_growth_model_types = ["ReactiveDimerization", "DimerCoalescence" ,"EBridgeModified", "IrreversibleDimerization"];

    n_sec = 60;
    sec_headers = ["N_agg", "d_m"];
    params_dict = {};
    for header in sec_headers:
        params_dict[header] = [f"{header}_{sec}" for sec in range(n_sec)];


    # Reading numerical results
    results_df = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                dict
            )
        )
    );

    for case in cases:
        for PAH_growth_model_type in PAH_growth_model_types:
            filename = f"results/{case}/KAUST/Sectional/{PAH_growth_model_type}/sim_results.csv";
            df = pd.read_csv(filename);
            results_df[case][PAH_growth_model_type] = df;

    # Reading experimental data
    PSD_df = {};
    for case in cases:
        PSD_df[case] = pd.read_csv(f"data/PSD_{case}.csv");
    N_agg_exp_df = pd.read_csv(f"data/N_agg.csv");


    z_sample = 1.4;

    # ------------------------------------------------------------------------------------------------------------------------
    # Particle size distribution
    # ------------------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 3, figsize=(figsize[0]*3+1, figsize[1]+0.5));


    for case_index, case in enumerate(cases):
        title = f"({string.ascii_lowercase[case_index]}) "+ make_title(case);
        # Numerical
        ax = axes[case_index];
    #     ax.set_title(title)
        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[case][PAH_growth_model_type];
            # Length
            z_array = df["z[m]"].to_numpy();
            i0 = np.where(z_array>z_sample)[0][0]-1;
            i1 = i0 + 1;
            w0 = (z_sample - z_array[i0])/(z_array[i1] - z_array[i0]);
            w1 = (z_array[i1] - z_sample)/(z_array[i1] - z_array[i0]); 
            # d_m
            params = [f"d_m_{i}" for i in range(60)];
            sub_df = df[params];
            d_m = ((sub_df.iloc[i0, :] * w0 + sub_df.iloc[i1, :] * w1)*1e9).to_numpy();
            x_data = (d_m[1:] + d_m[:-1])/2.0;
            # N_agg
            params = [f"N_agg_{i}" for i in range(60)];
            sub_df = df[params];
            #     =         mol/kg     * 1/mol  *     kg/m3  * 1e-6 m3/cc
            N_agg = ((sub_df.iloc[i0, :] * w0 + sub_df.iloc[i1, :] * w1) *   Av   *   gas_density / 1e6).to_numpy();
            y_data = ((N_agg[1:]+N_agg[:-1])/2.0)/np.diff(np.log(d_m));
            
            label = f"{labels[PAH_growth_model_type]}";
        
            ax.plot(x_data, y_data,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = label, linestyle = "solid");
        
        # Experimental
        x_data_exp = PSD_df[case]["d_m[nm]"];
        y_data_exp = PSD_df[case]["dN/dlogdm"];
        ax.plot(x_data_exp, y_data_exp, marker = "o", markersize = 12, markerfacecolor = "#FFFFFF",
            color="#000000", linewidth = linewidth, label = "Data", linestyle = "None");        

        
        # Y axis
        ax.set_yscale('log');
        ax.set_ylim([1e6, 3e11]);
        if case_index == 0:
            ax.set_ylabel("$dN_{agg}/d\\mathrm{log}(d_m)\\:[1/\\mathrm{cm}^3]$", {'fontsize' : axis_fontsize})
            

        # X axis
        ax.set_xscale('log');
        ax.set_xlabel("${d_m}$ [nm]", {'fontsize' : axis_fontsize})
        ax.set_xlim([1, 100]);
        
        title_pos = [
            {"x": 0.3, "y": 0.90},
            {"x": 0.3, "y": 0.90},
            {"x": 0.3, "y": 0.90},
        ]
        ax.text(title_pos[case_index]["x"], title_pos[case_index]["y"], title, transform=ax.transAxes, fontsize = title_fontsize);
        
        if case_index == 2:
            ax.legend(loc = (0.5, 0.2), ncol = 1, fontsize = legend_fontsize, framealpha=0, labelspacing=0.6, handletextpad=0.4);

    save_plot("PSD_diffQ", postprocess_dir = postprocess_dir, pad = 0.4);

    # ------------------------------------------------------------------------------------------------------------------------
    # N_agg
    # ------------------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 3, figsize=(figsize[0]*3+1, figsize[1]+0.9));

    for case_index, case in enumerate(cases):
        title = f"({string.ascii_lowercase[case_index]}) "+ make_title(case);
        # Numerical
        ax = axes[case_index];
        ax.set_title(title, pad = 13);
        
        param = "N_agg[#/cm3]";
        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[case][PAH_growth_model_type];
            # Length
            z_array = df["z[m]"].to_numpy();
            # N_agg
            f_array = df[param].to_numpy();

            label = labels[PAH_growth_model_type];        
            ax.plot(z_array, f_array,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = label, linestyle = "solid");
            
            
        
        # Experimental
        x_data_exp = 1.4;
        y_data_exp = N_agg_exp_df[case].to_numpy()[-1];
        ax.plot(x_data_exp, y_data_exp, marker = "s", markersize = 12, markerfacecolor = "#FFFFFF",
            color="#000000", linewidth = linewidth, label = "Data", linestyle = "None");        

        
        # Y axis
        ax.set_yscale('log');
        ax.set_ylim([1e7, 2e11]);
        if case_index == 0:
            ax.set_ylabel("$N_{agg} \\:[1/\\mathrm{cm}^3]$", {'fontsize' : axis_fontsize})
            
        

        # X axis
        ax.set_xlabel("z [m]", {'fontsize' : axis_fontsize});
        x0, x1 = 0, 1.42;
        ax.set_xlim([x0, x1]);
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(MultipleLocator(0.04));
        
        if case_index == 2:
            ax.legend(loc = "lower left", ncol = 1, fontsize = legend_fontsize, framealpha=0, labelspacing=0.4);
        
        
        param = "T[K]";
        # Numerical
        df = results_df[case][PAH_growth_model_type];
        # Length
        z_array = df["z[m]"].to_numpy();
        # N_agg
        f_array = df[param].to_numpy();
        
        
        z0 = z_array[f_array > 1673*0.9][0];
        z1 = z_array[f_array > 1673*0.9][-1];    
        rect = Rectangle((z0, 1), z1-z0, 1e14, linewidth=1, edgecolor='None', facecolor='#fc9403', alpha = 0.1);
        ax.add_patch(rect);
        
        ax.text((z0 + z1)/2 - 0.15, 8e10, "Hot zone", fontsize = title_fontsize, color = "#de0021");
    
    save_plot("N_agg", postprocess_dir = postprocess_dir, pad = 0.4);

    # ------------------------------------------------------------------------------------------------------------------------
    # Inception Flux
    # ------------------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 3, figsize=(figsize[0]*3+1, figsize[1]+0.9));
    for case_index, case in enumerate(cases):
        param = "inception_mass[mol/kg-s]"
        title = f"({string.ascii_lowercase[case_index]}) "+ make_title(case);
        # Numerical
        ax = axes[case_index];
        ax.set_title(title, pad = 13);
        for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
            df = results_df[case][PAH_growth_model_type];
            # Density
            rho_array = df["density[kg/m3]"].to_numpy();
            # Length
            z_array = df["z[m]"].to_numpy();
            # N_agg
            f_array = df[param].to_numpy() * Av * rho_array /1e6;

            label = labels[PAH_growth_model_type];        
            ax.plot(z_array, f_array,
                color=colors[PAH_growth_model_type_index], linewidth = linewidth, label = label, linestyle = "solid");
            
        # Y axis
        ax.set_yscale('log');
        ax.set_ylim([1e7, 4e12]);
        if case_index == 0:
            ax.set_ylabel("Soot inception flux, $I_{inc}$ [$1/\\mathrm{cm^3}-s$]", {'fontsize' : fontsize})
            

        # X axis
        ax.set_xlabel("z [m]", {'fontsize' : 22})
        x0, x1 = 0.0, 1.42;
        ax.set_xlim([x0, x1]);
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(MultipleLocator(0.04));
        
        legend_loc = [
            (0.36, 0.04),
            "center left",
            "center left"
        ]
        if case_index == 0:
            legend = ax.legend(loc = legend_loc[case_index], ncol = 1, fontsize = legend_fontsize,
                            framealpha=0, labelspacing=0.4, handlelength = 1.5);
        
        title_pos = [
            {"x": 0.3, "y": 0.90},
            {"x": 0.3, "y": 0.90},
            {"x": 0.3, "y": 0.90},
        ]        
        
        param = "T[K]";
        df = results_df[case][PAH_growth_model_type];
        # Length
        z_array = df["z[m]"].to_numpy();
        # N_agg
        f_array = df[param].to_numpy();
        # Hot zone
        z0 = z_array[f_array > 1673*0.9][0];
        z1 = z_array[f_array > 1673*0.9][-1];    
        rect = Rectangle((z0, 1), z1-z0, 1e14, linewidth=1, edgecolor='None', facecolor='#fc9403', alpha = 0.1);
        ax.add_patch(rect);
        
        ax.text((z0 + z1)/2 - 0.15, 1e12, "Hot zone", fontsize = title_fontsize, color = "#de0021");
    save_plot("inception", postprocess_dir = postprocess_dir, pad = 0.4);


    # ------------------------------------------------------------------------------------------------------------------------
    # Carbon mass from different processes
    # ------------------------------------------------------------------------------------------------------------------------
    fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize=(figsize[0]*2 - 0.5, figsize[1]*2 - 1.5));
    axes = axes.reshape(-1,)

    colors_set = ["#FE9D52", "#DDDDDD", "#2A6A99"];
    d_color = "#9c0000";
    d_color_graph = "#9c0000"

    case_index = 0
    case = cases[case_index];

    headers2 = ["d_p", "d_m"]
    headers2_dict = {"d_p": "d_p[nm]", "d_m": "d_m[nm]"}

    f_data_dict = {},
    for PAH_growth_model_type_index, PAH_growth_model_type in enumerate(PAH_growth_model_types):
        letter_index = PAH_growth_model_type_index;
        title = f"({string.ascii_lowercase[letter_index]}) " + labels[PAH_growth_model_type];
        ax = axes[PAH_growth_model_type_index];
        ax.text(0.4, 0.87, title, transform=ax.transAxes, fontsize = 23 , color = "#000000")

        header = "C_tot_int";
        df = results_df[case][PAH_growth_model_type];
        f_data_tot = df["C_tot_int[mol/kg]"];

        x_data = df["z[m]"]

        header = "C_tot_inception_int";
        f_data_inc = df["C_tot_inception_int[mol/kg]"];

        header = "C_tot_PAH_int";
        f_data_PAH = df["C_tot_PAH_int[mol/kg]"];

        header = "C_tot_growth_int";
        f_data_HACA = df["C_tot_growth_int[mol/kg]"];

        f_data_dict = {
            "Inception" : f_data_inc/f_data_tot,
            "PAH Adsorption" : f_data_PAH/f_data_tot,
            "HACA" : f_data_HACA/f_data_tot,
        }

        ax.stackplot(x_data, f_data_dict.values(), labels = f_data_dict.keys(), colors = colors_set,
                    alpha = 0.7);


        # Y axis
        if PAH_growth_model_type_index % 2 == 0:
            ax.set_ylabel("$\\frac{\\mathrm{Soot\\;carbon\\;mass\\;from\\;each\\;source}}{\\mathrm{Soot\\;carbon\\;mass}}$",
                        {'fontsize' : axis_fontsize+2});
        ax.set_ylim([0, 1]);
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.04))

        # X axis
        if PAH_growth_model_type_index > 1:
            ax.set_xlabel("z [m]", {'fontsize' : 22});
        ax.tick_params(axis='x', pad = 10)
        x0, x1 = 0.4, 1.4;
        ax.set_xlim([x0, x1]);
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(MultipleLocator(0.04));


        legend_pos = [
            (0.027, 0.4),
            (0.027, 0.4)
        ];


        if PAH_growth_model_type_index == 3:
            legend = ax.legend(loc = (0.36, 0.5), ncol = 1, fontsize = 20, columnspacing = 0.1,
                            framealpha = 1, handletextpad = 0.3, handlelength = 0.8);


        ax2 = ax.twinx();

        # Length
        z_array = df["z[m]"].to_numpy();
        for header_index, header in enumerate(headers2):
            param = headers2_dict[header];
            f_array = df[param].to_numpy();
            label = labels[header]
            ax2.plot(z_array, f_array,
                color = d_color_graph, linewidth = linewidth, label = label,
                    linestyle = linestyles[header_index]);
        if PAH_growth_model_type_index % 2 == 1:
            ax2.set_ylabel("Diameter [nm]")
        ax2.set_ylim([1, 30]);
        ax2.yaxis.set_major_locator(MultipleLocator(5))
        ax2.yaxis.set_minor_locator(MultipleLocator(1));


        ax2.spines['right'].set_color(d_color);
        ax2.tick_params(axis='y', colors=d_color)
        ax2.yaxis.label.set_color(d_color);

        if PAH_growth_model_type_index == 3:
            legend = ax2.legend(loc = (0.5, 0.2), ncol = 1, fontsize = 20, columnspacing = 0.1,
                            framealpha = 0.0, handletextpad = 1, handlelength = 2);
            for text in legend.get_texts():
                text.set_color(d_color_graph)

    save_plot("C_tot_distmap", postprocess_dir = postprocess_dir, pad = 0.4)

if __name__ == "__main__":
    plot();