import os
import shutil
import time
import numpy as np

from wisdem import run_wisdem
from wisdem.inputs import load_yaml, write_yaml

# User Inputs for driver

# Should we go through and run the design optimizations?
run_optimization = True
overwrite_geometry = True

# Should we just evaluate the current design and generate WISDEM outputs? (always run if optimizing)
run_evaluation = True

# Plot results
plot_results = False

# Set up paths
CWDDIR = os.path.dirname(os.path.realpath(__file__))
fanalysis = os.path.join(CWDDIR, 'analysis_options.yaml')
fanalysis_opt = os.path.join(CWDDIR, 'analysis_options_monopile.yaml')
fmodeling = os.path.join(CWDDIR, 'modeling_options_monopile.yaml')
fanalysis_tmp = 'analysis_options_temp.yaml'
analysis_opt_yaml = load_yaml(fanalysis_opt)

rating = 15
depths = [20] #, 30, 40, 50, 60]
# Set the maximum diameter [in m] for the optimizations.  Can be constant or refine by rating-depth combo
max_diam = 12 * np.ones(len(depths))

# Set the first natural frequency [in Hz] of the monopile-tower structure (with the nacelle+rotor on top)
freq_lower = 0.22 * np.ones(len(depths))
freq_upper = 0.4 * np.ones(len(depths))
# ---------------------------------------------------------------------------------------------

if True:
    fgeometry = f'modified_{rating}.yaml'
else:
    fgeometry = f'modified_{rating}_gb.yaml'

# The shortened WISDEM simulation that only runs the tower-monopile modules for quicker optimization
fmodeling_opt  = os.path.join(CWDDIR, f'{rating}mw', 'modeling_options_monopile_noRNA.yaml')

os.chdir(os.path.join(CWDDIR, f'{rating}mw'))
for di, d in enumerate(depths):

    os.chdir(f'{d}m')
    if run_optimization:
    # Write out customized analysis options
        analysis_opt_yaml['design_variables']['tower']['outer_diameter']['upper_bound'] = float(max_diam[di])
        analysis_opt_yaml['design_variables']['monopile']['outer_diameter']['upper_bound'] = float(max_diam[di])
        analysis_opt_yaml['constraints']['tower']['frequency_1']['lower_bound'] = float(freq_lower[di])
        analysis_opt_yaml['constraints']['tower']['frequency_1']['upper_bound'] = float(freq_upper[di])
        analysis_opt_yaml['constraints']['monopile']['frequency_1']['lower_bound'] = float(freq_lower[di])
        analysis_opt_yaml['constraints']['monopile']['frequency_1']['upper_bound'] = float(freq_upper[di])
        write_yaml(analysis_opt_yaml, fanalysis_tmp)

        # Run WISDEM optimization of tower and monopile
        t = time.time()
        #wt_opt, _, _ = run_wisdem(fgeometry, fmodeling_opt, fanalysis_tmp)

        # Read output
        fopt_path = os.path.join('outputs', 'monotow_opt_output.yaml')
        if not os.path.exists(fopt_path): continue

        # Overwrite orignal file
        if overwrite_geometry:
            shutil.move(fopt_path, fgeometry)

        print(f"{r}mw {d}m Full Opt Elapsed Time: {time.time() - t:.2f}")

    if run_evaluation or run_optimization:

        wt_run, _, _ = run_wisdem(fgeometry, fmodeling, fanalysis)

    os.chdir('..')