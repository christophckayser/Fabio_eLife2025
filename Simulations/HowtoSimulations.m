% Howto for the simulations of effects


% ------------------------------------------------------------
% test analysis for single trial data
% ------------------------------------------------------------
Trials_run
% calls ... to work
Trials_simulate

% calls ... to generate data
local_sim_gendata
% calls ... for analysis
local_fitmodelsTrial


% ------------------------------------------------------------
% test analysis for binned data
% does both spectral approach and linear model based approach
% ------------------------------------------------------------
Trials_run
% calls ... to work
Binned_simulate

% calls ... to generate data
local_sim_gendata
% calls ... for analysis
local_fitmodelsGLM
local_computespectra


% ------------------------------------------------------------
% visualize results = Figure 2
% ------------------------------------------------------------
CK15_MakeFigure2