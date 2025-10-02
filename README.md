### Code for 'The membrane transition strongly enhances biopolymer condensation through prewetting' ###
## Simulation

### Compile the code with the command 
  'python setup_sim.py build_ext --inplace' 
### Run the Simulations with
  'python run_prewetting_sim.py [c_tether] [J_bulk] [chem_potential_1] [chem_potential_2] [membrane_composition] [Tmembrane/Tc] [Idx] ' 
  Command line parmaaters are
  # ctether - tether concentration
  # J_bulk - contact coupling between molecules
  # chemical potential 1 - chemical potential of species 1
  # chemical potential 2 - chemical potential of species 2
  # membrane composition - ratio of up spins / total spins of membrane (0 = single component)
  # Tmembrane / Tc - Interactions between membrane spins
  # idx - indicator variable for replicates

## Generate commands with
'python write_jobs_figs.py run_prewetting_sim.py' 
This will generate *.sh files with example commands, that were used to generate the data in the figures. 

### Plotting
Plot and analyze configuration files with code in this directory
# generate_snapshots_prewetting.py - generates sample configurations. Inputs:  a file that lists configuration file name formatted as $config_bulk\t$config_surface\n
# plot_prewetting_exp_dense_prof - genrates a lateral 'density profile' of bulk. Inputs: a file that lists configuration file name formatted as $config_bulk\t$config_surface\n
# plot_prewetting_exp_sim_data_coex.py - generats *ads_data.txt file of the adsorption, phase area fractions. Inputs $config_bulk 
# plot_prewetting_exp_width - generates a file of fit widths of preweti domains. Inputs: a file that lists configuration file name formatted as $config_bulk\t$config_surface\n


