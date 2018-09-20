# ACC2019
Matlab code for Paper submmited to the American Control Conference 2019.
Required: Mosek: SDP package and Yalmip
data.txt: File with experimental data, native
crossover_data.mat: experimental data in Matlab format
crossover_lmi_lme_design.m: obtain observer gains --> obs_gains.mat
crossover_lmi_lme_polytopic_design.m: obtain observer gains via polytopic solution of the LMI-LME problem--> obs_gains_poly.mat
crossover_sys.m: simulation for linear crossover dynamic
crossover_sys_data.m: simulation for linear crossover dynamic and comparison with experimental data
crossover_obs_sys.m: Simulation Adaptive Observer and system with linear crossover dynamic
crossover_obs_data.m: Simulation Adaptive Observer, system with linear crossover dynamic and experimetal data
