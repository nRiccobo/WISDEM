# Needed libraries for this script
import os
import shutil
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

import matplotlib.pyplot as plt
from wisdem import run_wisdem
from helper_functions import load_yaml, write_yaml, init_fig_axis, fig_format, save

import time

# ---------------------------------------------------------------------------------------------
### USER OPTIONS FOR SCRIPT ###

# Should we go through and run the design optimizations?
run_optimization = True #True
overwrite_geometry = True

reset_geometry = False

# Should we just evaluate the current design and generate WISDEM outputs? (always run if optimizing)
run_evaluation = True

# Is the goal to evaluate support fatigue with DLC 1.2
run_fatigue_openfast = False

# Should we generate output plots
make_plots = False

# Gearbox or Direct-drive
gear_box = False

# Which machine ratings should we run [in MW] (choices are 15, 20, 22, 25)?
ratings = [15] # , 20, 22] # [15, 20, 22, 25]

# Which depths should we run [in m] (choices are 20, 30, 40, 50, 60)?
depths  = [20, 30, 40, 50, 60] # [20, 30, 40, 50, 60]

# Set the maximum diameter [in m] for the optimizations.  Can be constant or refine by rating-depth combo
max_diam = 13 * np.ones( (len(ratings), len(depths)) )
tow_diam = 10 * np.ones( (len(ratings), len(depths)) )
#if len(ratings) > 1:
#    max_diam[1,:] = 11. # 20 m
#    max_diam[2,:] = 11. # 22 m
    #max_diam[3,:] = 12. # 25 m

# Set the first natural frequency [in Hz] of the monopile-tower structure (with the nacelle+rotor on top)
freq_lower = 0.13 * np.ones( (len(ratings), len(depths)) )
freq_upper = 0.24 * np.ones( (len(ratings), len(depths)) )
# ---------------------------------------------------------------------------------------------

# Load WEIS only if that is what is being run, otherwise it brings in many new dependencies and baggage
if run_fatigue_openfast:
    from weis.glue_code.runWEIS import run_weis

# Set paths to input files starting with the directory for this file as the 'root'
mydir         = os.path.dirname(os.path.realpath(__file__))  # get path to this file
fanalysis     = os.path.join(mydir, 'analysis_options.yaml')
fanalysis_opt = os.path.join(mydir, 'analysis_options_monopile.yaml')
fmodeling     = os.path.join(mydir, 'modeling_options_monopile.yaml')
fanalysis_tmp = 'analysis_options_temp.yaml'
ftmp          = 'optimized.yaml'

# Load in optimization analysis options for edits later
analysis_opt_yaml = load_yaml(fanalysis_opt)

# Loop over all ratings, and descend into folder
for ri, r in enumerate(ratings):

    # Declare WISDEM and WEIS input files that vary by machine rating:
    # The turbine definition file
    if gear_box:
        fgeometry = f'modified_{r}_gen.yaml'

    else:
        fgeometry = f'modified_{r}.yaml'

    # The shortened WISDEM simulation that only runs the tower-monopile modules for quicker optimization
    fmodeling_opt  = os.path.join(mydir, f'{r}mw', 'modeling_options_monopile_noRNA.yaml')

    os.chdir(os.path.join(mydir, f'{r}mw'))

    for di, d in enumerate(depths):

        os.chdir(f'{d}m')

        # Run optimization
        if run_optimization:

            # reset tower and monopile geometry
            if reset_geometry:
                write_yaml(fgeometry, fanalysis_tmp)
                fgeometry_tmp = load_yaml(fgeometry)

                fgeometry_tmp['components']['monopile']['outer_shape_bem']['outer_diameter']['values'] = [9.5]*7
                #fgeometry_tmp['components']['monopile']['internal_structure_2d_fem']['layers']['thickness']['values'] = [0.09]*6

            else:
                fgeometry_tmp = fgeometry

            # Write out customized analysis options
            analysis_opt_yaml['design_variables']['tower']['outer_diameter']['upper_bound'] = float(tow_diam[ri, di])
            analysis_opt_yaml['design_variables']['monopile']['outer_diameter']['upper_bound'] = float(max_diam[ri, di])
            analysis_opt_yaml['constraints']['tower']['frequency_1']['lower_bound'] = float(freq_lower[ri, di])
            analysis_opt_yaml['constraints']['tower']['frequency_1']['upper_bound'] = float(freq_upper[ri, di])
            analysis_opt_yaml['constraints']['monopile']['frequency_1']['lower_bound'] = float(freq_lower[ri, di])
            analysis_opt_yaml['constraints']['monopile']['frequency_1']['upper_bound'] = float(freq_upper[ri, di])
            write_yaml(analysis_opt_yaml, fanalysis_tmp)

            # Run WISDEM optimization of tower and monopile
            t = time.time()
            wt_opt, _, _ = run_wisdem(fgeometry_tmp, fmodeling_opt, fanalysis_tmp)

            # Read output
            fopt_path = os.path.join('outputs_mono', 'monotow_iea15_output.yaml')
            if not os.path.exists(fopt_path): continue

            # Overwrite orignal file
            if overwrite_geometry:
                shutil.move(fopt_path, fgeometry)

            print(f"{r}mw {d}m Full Opt Elapsed Time: {time.time() - t:.2f}")

        # Run full WISDEM
        if run_evaluation or run_optimization:

            t2 = time.time()
            wt_run, _, _ = run_wisdem(fgeometry, fmodeling, fanalysis)

            print(f"{r}mw {d}m Eval Elapsed Time: {time.time() - t2:.2f}")

        # Run WEIS to calculate fatigue results
        if run_fatigue_openfast:
            wt_run, _, _ = run_weis(fgeometry, fmodeling, fanalysis)

        os.chdir('..')

    os.chdir('..')

if make_plots:
    # Data containers to hold outputs from test matrix
    m_tow = np.zeros((len(ratings), len(depths)))
    c_tow = np.zeros(m_tow.shape)
    m_mono = np.zeros(m_tow.shape)
    c_mono = np.zeros(m_tow.shape)
    m_trans = np.zeros(m_tow.shape)
    c_trans = np.zeros(m_tow.shape)
    m_struct = np.zeros(m_tow.shape)
    rating_mat = np.zeros(m_tow.shape)
    depth_mat = np.zeros(m_tow.shape)
    m_nacelle = np.zeros(m_tow.shape)
    c_turbine = np.zeros(m_tow.shape)

    # When there is an array of data at every point, that could change size, use nested dictionaries instead of 3-D matrices
    zpts = {}
    diam = {}
    thick = {}

    # Loop over all ratings, and descend into folder
    for ri, r in enumerate(ratings):
        os.chdir(f'{r}mw')

        # Initialize nested dictionaries
        zpts[r] = {}
        diam[r] = {}
        thick[r] = {}

        # Loop over all depths, and descend into folder
        for di, d in enumerate(depths):
            os.chdir(f'{d}m')

            # Store ratings and depths for easy plotting
            rating_mat[ri,di] = r
            depth_mat[ri,di] = d

            # Read archive of outputs from WISDEM run
            fout_path = os.path.join('outputs', 'equinor_turb_output.npz')
            if not os.path.exists(fout_path): continue
            out_archive = np.load( fout_path )

            # Store outputs into matrices.  List of available outputs can be found in the WISDEM documentation or in the csv- or xlsx-files
            m_tow[ri,di] = out_archive['towerse.tower_mass_kg']
            c_tow[ri,di] = out_archive['towerse.tower_cost_USD']
            m_mono[ri,di] = out_archive['fixedse.monopile_mass_kg']
            c_mono[ri,di] = out_archive['fixedse.monopile_cost_USD']
            m_struct[ri,di] = out_archive['fixedse.structural_mass_kg']
            m_nacelle[ri,di] = out_archive['orbit.nacelle_mass_t']
            c_turbine[ri,di] = out_archive['orbit.turbine_capex_USD/kW']

            # Store the arrays for diameter, thickness, and z-points in the nested dictionary
            zpts[r][d] = np.r_[out_archive['monopile.ref_axis_m'][:,2], out_archive['tower.ref_axis_m'][:,2]]
            diam[r][d] = np.r_[out_archive['monopile.diameter_m'], out_archive['tower.diameter_m']]
            thick[r][d] = np.r_[out_archive['monopile.layer_thickness_m'][0,:], out_archive['tower.layer_thickness_m'][0,:]]

            os.chdir('..')
        os.chdir('..')

    # Create figures of design point trends for tower mass, monopile mass, combined mass
    xstr = 'Depth [m]'
    legstr = ['15 MW','20 MW','22 MW','25 MW']
    fig, ax = init_fig_axis()
    ax.plot(depths, 1e-3*m_tow.T, '-o')
    ax.grid()
    ax.set_xlabel(xstr)
    ax.set_ylabel('Tower mass [t]')
    ax.legend(legstr)
    fig_format(ax)
    save(fig, 'mono_cases-tower_mass')

    fig, ax = init_fig_axis()
    ax.plot(depths, 1e-3*m_mono.T, '-o')
    ax.grid()
    ax.set_xlabel(xstr)
    ax.set_ylabel('Monopile mass [t]')
    ax.legend(legstr)
    fig_format(ax)
    save(fig, 'mono_cases-monopile_mass')

    fig, ax = init_fig_axis()
    ax.plot(depths, 1e-3*m_struct.T, '-o')
    ax.grid()
    ax.set_xlabel(xstr)
    ax.set_ylabel('Support structure mass [t]')
    ax.legend(legstr)
    fig_format(ax)
    save(fig, 'mono_cases-structural_mass')

    # Tower geometry profiles by depth
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    brown = np.array([150.0, 75.0, 0.0]) / 256.0
    for di, d in enumerate(depths):
        ftow = plt.figure(figsize=(11, 4))
        ax1 = ftow.add_subplot(121)
        for ri, r in enumerate(ratings):
            ax1.plot(diam[r][d], zpts[r][d], '-', color=colors[ri], label=f'{r} MW')
        vx = ax1.get_xlim()
        zs = zpts[r][d]
        if zs.min() < -5.0:
            h_trans = out_archive['monopile.ref_axis_m'][-1,2]
            ax1.plot(vx, np.zeros(2), color="b", linestyle="--")
            ax1.plot(vx, -d * np.ones(2), color=brown, linestyle="--")
            ax1.plot(vx, h_trans * np.ones(2), color="k", linestyle="--")
            ax1.text(vx[0] + 0.02 * np.diff(vx), 2, "Water line", color="b", fontsize=12)
            ax1.text(vx[0] + 0.02 * np.diff(vx), -d + 2, "Mud line", color=brown, fontsize=12)
            ax1.text(vx[0] + 0.02 * np.diff(vx), h_trans + 2, "Tower transition", color="k", fontsize=12)
        ax1.set_xlim(vx)
        plt.xlabel("Outer Diameter [m]")
        plt.ylabel("Tower Height [m]")
        plt.grid(color=[0.8, 0.8, 0.8], linestyle="--")

        ax2 = ftow.add_subplot(122)
        for ri, r in enumerate(ratings):
            ax2.step(1e3 * thick[r][d], zpts[r][d], '-', color=colors[ri], label=f'{r} MW', where="post")
        vx = ax2.get_xlim()
        if zs.min() < -5.0:
            ax2.plot(vx, np.zeros(2), color="b", linestyle="--")
            ax2.plot(vx, -d * np.ones(2), color=brown, linestyle="--")
            ax2.plot(vx, h_trans * np.ones(2), color="k", linestyle="--")
        ax2.set_xlim(vx)
        plt.xlabel("Wall Thickness [mm]")
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.grid(color=[0.8, 0.8, 0.8], linestyle="--")
        plt.subplots_adjust(bottom=0.15, left=0.15)
        ftow.subplots_adjust(hspace=0.02, wspace=0.02, bottom=0.15, left=0.15)
        save(ftow, f'tower-monopile_depth{d}') #, pad_inches=0.1, bbox_inches="tight")

