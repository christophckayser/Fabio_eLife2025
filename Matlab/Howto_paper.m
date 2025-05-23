

% Analysis pipeline to reproduce the analyses and figures

% set directories. needs to be adapted accordingly
defdirsCK15 

% /func  contains functions called by the analysis routines. They are
% mentioned below where needed

% /simulations contains the code for running the stimulations

%------------------------------------------------------
%  Figure 3 - Visualize group-level data
%------------------------------------------------------
CK15_MakeFigure3

% relies on ... to compute the SOA binned behavior metrics
compute_binned_behavior 


% --------------------------------------------------------
% implement AR-based analysis of spectra
% --------------------------------------------------------
Analyze_CK15_Exp1234_spectra
Analyze_CK15_Exp1234_spectra_display

% relies on ... to compute the SOA binned behavior metrics
compute_binned_behavior 

% relies on ... to compute participant-wise spectra and AR surrogates
local_computespectra
% requires arfit package!


% --------------------------------------------------------
% implement regression analysis of binned data 
%------------------------------------------------------
Analyze_CK15_Exp1234_bins

% relies on ... to fit GLM models to participant wise data
local_fitmodelsGLM
% relies on 
fast_regress
ck_decoding_logist


%------------------------------------------------------
%  Figure 4 - Visualize group-level statistical results
%------------------------------------------------------
CK15_MakeFigure4


CK15_MakeFigure5

% --------------------------------------------------------
% implement simulations of participant samples to compute prevalences
% --------------------------------------------------------
Analyze_CK15_Exp1234_spectra_MC
Analyze_CK15_Exp1234_bins_MC
Analyze_CK15_Exp123_singletrial_MC

%------------------------------------------------------
%  Figure 5/6 - Prevalence results
%------------------------------------------------------
CK15_MakeFigure6and7


% ------------------------------------------------------------
% analysis of eye tracking data in Experiment 4
% ------------------------------------------------------------
% analysis of Eye movements. needs to be run before the below
Analyze_CK15_Exp4_collecteyedata

%------------------------------------------------------
%  Figure 9 eye movement data
%------------------------------------------------------
CK15_MakeFigure9
% not used in the paper anymore


% ------------------------------------------------------------
% analysis of behavioral data Experiment 4
% ------------------------------------------------------------
% data split by conditions
Analyze_CK15_Exp4_conditions

% data split by eye movement properties
Analyze_CK15_Exp4_eyesplits
% requires result from Analyze_CK15_Exp4_eyedata.m

% relies on ... to fit GLM models to participant wise data
local_fitmodelsGLM
% relies on 
fast_regress
ck_decoding_logist


% for these we also compute the prevalences
Analyze_CK15_Exp4_conditions_MC
Analyze_CK15_Exp4_eyesplits_MC


%------------------------------------------------------
%  Figure 8 results exp 4
%------------------------------------------------------
CK15_MakeFigure8



% --------------------------------------------------------
% method simulations
% F:\CKDATA\Projects\projects\Hearing\CK15\matlab\simulations
Display_SimResult
