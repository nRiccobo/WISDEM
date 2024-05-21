# Needed libraries for this script
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from wisdem import run_wisdem
from helper_functions import load_yaml, write_yaml, init_fig_axis, fig_format, save

# ---------------------------------------------------------------------------------------------
### USER OPTIONS FOR SCRIPT ###

# Should we go through and run the design optimizations?
run_optimization = False #True
overwrite_geometry = False

# Should we just evaluate the current design and generate WISDEM outputs? (always run if optimizing)
run_evaluation = False

# Should we generate output plots
make_plots = True

# Which machine ratings should we run [in MW] (choices are 15, 20, 22, 25)?
ratings = [15, 20, 22, 25]

# Which depths should we run [in m] (choices are 20, 30, 40, 50, 60)?
depths  = [20, 30, 40, 50, 60]

# Number of legs to consider in jacket design (choices are 3, 4)?
legs = [3] #, 4]

# Set the maximum diameter [in m] for the optimizations.  Can be constant or refine by rating-depth combo
max_head_radius = 10. * np.ones( (len(ratings), len(depths)) )
max_head_radius[1,:] = 11. # 20 m
max_head_radius[2,:] = 11. # 22 m
max_head_radius[3,:] = 12. # 25 m

# Set the first natural frequency [in Hz] of the tower (with the nacelle+rotor on top)
freq_lower = 0.13 * np.ones( (len(ratings), len(depths)) )
freq_upper = 0.4 * np.ones( (len(ratings), len(depths)) )
# ---------------------------------------------------------------------------------------------

# Set paths to input files starting with the directory for this file as the 'root'
mydir         = os.path.dirname(os.path.realpath(__file__))  # get path to this file
fanalysis     = os.path.join(mydir, 'analysis_options.yaml')
fanalysis_opt = os.path.join(mydir, 'analysis_options_jacket.yaml')
fmodeling     = os.path.join(mydir, 'modeling_options_jacket.yaml')
fanalysis_tmp = 'analysis_options_temp.yaml'
ftmp          = 'optimized.yaml'

# Load in optimization analysis options for edits later
analysis_opt_yaml = load_yaml(fanalysis_opt)

# Loop over all number of leg options, and descend into folder
for lgi, lg in enumerate(legs):
    os.chdir(f'{lg}leg')

    # Loop over all ratings, and descend into folder
    for ri, r in enumerate(ratings):
        # Declare WISDEM and WEIS input files that vary by machine rating:
        # The turbine definition file
        fgeometry = f'modified_{r}.yaml'

        # The shortened WISDEM simulation that only runs the tower-jacket modules for quicker optimization
        fmodeling_opt  = os.path.join(mydir, f'{lg}leg', f'{r}mw', 'modeling_options_jacket_noRNA.yaml')

        os.chdir(f'{r}mw')

        for di, d in enumerate(depths):
            os.chdir(f'{d}m')

            # Run optimization
            if run_optimization:
                # Write out customized analysis options
                analysis_opt_yaml['design_variables']['tower']['outer_diameter']['upper_bound'] = float(max_head_radius[ri, di])
                analysis_opt_yaml['design_variables']['jacket']['r_head']['upper_bound'] = float(max_head_radius[ri, di])
                analysis_opt_yaml['constraints']['tower']['frequency_1']['lower_bound'] = float(freq_lower[ri, di])
                analysis_opt_yaml['constraints']['tower']['frequency_1']['upper_bound'] = float(freq_upper[ri, di])
                write_yaml(analysis_opt_yaml, fanalysis_tmp)

                # Run WISDEM optimization of tower and jacket
                wt_opt, _, _ = run_wisdem(fgeometry, fmodeling_opt, fanalysis_tmp)

                # Read output
                fopt_path = os.path.join('outputs', 'jacktow_output.yaml')
                if not os.path.exists(fopt_path): continue

                # Overwrite orignal file
                if overwrite_geometry:
                    shutil.move(fopt_path, fgeometry)

            # Run full WISDEM
            if run_evaluation or run_optimization:
                wt_run, _, _ = run_wisdem(fgeometry, fmodeling, fanalysis)

            os.chdir('..')
        os.chdir('..')
    os.chdir('..')

if make_plots:
    # Data containers to hold outputs from test matrix
    m_tow = np.zeros((len(legs), len(ratings), len(depths)))
    c_tow = np.zeros(m_tow.shape)
    m_jack = np.zeros(m_tow.shape)
    c_jack = np.zeros(m_tow.shape)
    m_trans = np.zeros(m_tow.shape)
    c_trans = np.zeros(m_tow.shape)
    m_struct = np.zeros(m_tow.shape)
    rating_mat = np.zeros(m_tow.shape)
    depth_mat = np.zeros(m_tow.shape)
    leg_mat = np.zeros(m_tow.shape)
    rfoot = np.zeros(m_tow.shape)
    m_nacelle = np.zeros(m_tow.shape)
    c_turbine = np.zeros(m_tow.shape)

    # When there is an array of data at every point, that could change size, use nested dictionaries instead of 3-D matrices
    zpts = {}
    diam = {}
    thick = {}

    # Loop over all number of leg options, and descend into folder
    for lgi, lg in enumerate(legs):
        os.chdir(f'{lg}leg')
        
        zpts[lg] = {}
        diam[lg] = {}
        thick[lg] = {}

        # Loop over all ratings, and descend into folder
        for ri, r in enumerate(ratings):
            os.chdir(f'{r}mw')

            # Initialize nested dictionaries
            zpts[lg][r] = {}
            diam[lg][r] = {}
            thick[lg][r] = {}

            # Loop over all depths, and descend into folder
            for di, d in enumerate(depths):
                os.chdir(f'{d}m')

                # Store ratings and depths for easy plotting
                rating_mat[lgi,ri,di] = r
                depth_mat[lgi,ri,di] = d
                leg_mat[lgi,ri,di] = lg

                # Read archive of outputs from WISDEM run
                fout_path = os.path.join('outputs', 'equinor_turb_output.npz')
                if not os.path.exists(fout_path): continue
                out_archive = np.load( fout_path )

                # Store outputs into matrices.  List of available outputs can be found in the WISDEM documentation or in the csv- or xlsx-files
                m_tow[lgi,ri,di] = out_archive['towerse.tower_mass_kg']
                c_tow[lgi,ri,di] = out_archive['towerse.tower_cost_USD']
                m_jack[lgi,ri,di] = out_archive['fixedse.jacket_mass_kg']
                c_jack[lgi,ri,di] = out_archive['fixedse.jacket_cost_USD']
                m_struct[lgi,ri,di] = out_archive['fixedse.structural_mass_kg']
                m_nacelle[lgi,ri,di] = out_archive['orbit.nacelle_mass_t']
                c_turbine[lgi,ri,di] = out_archive['orbit.turbine_capex_USD/kW']

                # Store the arrays for diameter, thickness, and z-points in the nested dictionary
                zpts[lg][r][d] = np.r_[-d, out_archive['tower.ref_axis_m'][0,2], out_archive['tower.ref_axis_m'][:,2]]
                r_head = float(out_archive['jacket.r_head_m'])
                rfoot[lgi,ri,di] = r_foot = r_head*float(out_archive['jacket.foot_head_ratio'])
                diam[lg][r][d] = np.r_[2*r_foot, 2*r_head, out_archive['tower.diameter_m']]
                # Unsure how to plot thicknesses for the jacket
                #thick[lg][r][d] = np.r_[out_archive['jacket.layer_thickness_m'][0,:], out_archive['tower.layer_thickness_m'][0,:]]

                os.chdir('..')
            os.chdir('..')
        os.chdir('..')

    # Create figures of design point trends for tower mass, jacket mass, combined mass
    xstr = 'Depth [m]'
    legstr = ['15 MW','20 MW','22 MW','25 MW']
    fig, ax = init_fig_axis()
    ax.plot(depths, 1e-3*m_tow[0,:,:].T, '-o')
    if len(legs) > 1:
        ax.set_prop_cycle(None)
        ax.plot(depths, 1e-3*m_tow[1,:,:].T, '--o')
    ax.grid()
    ax.set_xlabel(xstr)
    ax.set_ylabel('Tower mass [t]')
    ax.legend(legstr)
    fig_format(ax)
    save(fig, 'jack_cases-tower_mass')

    fig, ax = init_fig_axis()
    ax.plot(depths, 1e-3*m_jack[0,:,:].T, '-o')
    if len(legs) > 1:
        ax.set_prop_cycle(None)
        ax.plot(depths, 1e-3*m_jack[1,:,:].T, '--o')
    ax.grid()
    ax.set_xlabel(xstr)
    ax.set_ylabel('Jacket mass [t]')
    ax.legend(legstr)
    fig_format(ax)
    save(fig, 'jack_cases-jacket_mass')

    fig, ax = init_fig_axis()
    ax.plot(depths, 1e-3*m_struct[0,:,:].T, '-o')
    if len(legs) > 1:
        ax.set_prop_cycle(None)
        ax.plot(depths, 1e-3*m_struct[1,:,:].T, '--o')
    ax.grid()
    ax.set_xlabel(xstr)
    ax.set_ylabel('Support structure mass [t]')
    ax.legend(legstr)
    fig_format(ax)
    save(fig, 'jack_cases-structural_mass')
    
    # Tower geometry profiles by depth
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    brown = np.array([150.0, 75.0, 0.0]) / 256.0
    for di, d in enumerate(depths):
        fig, ax = init_fig_axis()
        for ri, r in enumerate(ratings):
            if 3 in legs:
                ax.plot(diam[3][r][d], zpts[3][r][d], '-', color=colors[ri], label=f'{r} MW')
            if 4 in legs:
                ax.set_prop_cycle(None)
                ax.plot(diam[4][r][d], zpts[4][r][d], '--', color=colors[ri], label=f'{r} MW')
        vx = ax.get_xlim()
        zs = zpts[3][r][d]
        if zs.min() < -5.0:
            h_trans = out_archive['tower.ref_axis_m'][0,2]
            ax.plot(vx, np.zeros(2), color="b", linestyle="--")
            ax.plot(vx, -d * np.ones(2), color=brown, linestyle="--")
            ax.plot(vx, h_trans * np.ones(2), color="k", linestyle="--")
            ax.text(vx[0] + 0.02 * np.diff(vx), 2, "Water line", color="b", fontsize=12)
            ax.text(vx[0] + 0.02 * np.diff(vx), -d + 2, "Mud line", color=brown, fontsize=12)
            ax.text(vx[0] + 0.02 * np.diff(vx), h_trans + 2, "Tower transition", color="k", fontsize=12)
        ax.set_xlim(vx)
        plt.xlabel("Outer Diameter [m]")
        plt.ylabel("Tower Height [m]")
        plt.grid(color=[0.8, 0.8, 0.8], linestyle="--")
        fig_format(ax)
        save(fig, f'tower-jacket_depth{d}') #, pad_inches=0.1, bbox_inches="tight")

