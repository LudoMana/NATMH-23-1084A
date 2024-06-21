# README for NATMH-23-1084A study codes

# Empirical analysis:

run Empirical_Analysis.m in matlab to:
- load and organise data into subgroups using STAGE_vec.mat (128 healthy controls;37 stage 2 patients;31 stage3 remitting-relapsing patients; 20 stage3 non-remitting patients)
- compute global and local measures of pairwise functional connectivity and functional strength for each condition
- save empirical data

run Empirical_comparison_FC.m in matlab to:
- compare functional connectivity measures between conditions and compute statistics
- plot results
- save empirical statistics

# Model based analysis:

run LUDO_heter_SC_cnt_model_psycho_CNT_cluster.m in matlab on a cluster (script parallelised for iterations) to:
- set model details and fix parameters
- compute optimal vector of Bifurcation Parameters (a) for each value of global coupling (G) for healty controls
- simulate data and compare to empirical data of healthy controls (HC) to compute optimal fit
- save simulated data

run the equivalent versions of the above model script adapted to compute optimal fit for Early psychosis patients:
- LUDO_heter_SC_ep_model_psycho_EP3rem_rel_cluster.m for stage 3 remitting-relapsing patients (EP3R);
- LUDO_heter_SC_ep_model_psycho_EP3Nrem_cluster.m for stage 3 non remitting patients (EP3NR);
- LUDO_heter_SC_ep_model_psycho_EP2_cluster.m for stage 2 patients (EP2)

run LUDO_homo_SC_cnt_model_psycho_CNT_cluster in matlab on a cluster (script parallelised for iterations) to:
- compute the homogeneous version of the model for healtly controls (Bifurcation paramenter fixed to bifurcation point a=0 for all brain nodes)

run visualize_model_results in matlab to:
- aggregate results for all iterations for each condition
- plot model results for healthy controls
- compare optimal values of Global coumpling (G) and Bifurcation Parameters (a) between condition and compute statistics
- plot comparisons between optimal parameters in patients (EP3R, EP3NR, EP2) vs healty controls (HC)

run TOY_MODEL.m to:
- set model details and fix parameters for the simplified network model (11 nodes, 2 communities, 1 hub, 1 disconnected node)
- run simulations at rest at different values of Bifurcation Parameter associated to the hub and/or the disconnected node
- run simulations for increasing values incoming stimulation at different values of Bifurcation Parameter associated to the hub and/or the disconnected node
- compute pairwise and global values of network functional connectivity for each simulation
- plot results
