
% simulate different scenarious of rhythmicity in behavior
% generate data as sum of 5 components with fixed and random 
% beta's across participants. 
%
% testing binned analysis
% spectral analysis

clear;

%----------------------------------------------

ARG.DT = 60;  
ARG.BINS = [0:ARG.DT :1200]./1000; 
ARG.nbins = length(ARG.BINS)-1;

ARG.flist = [1.2:0.1:4, 4.2:0.2:8]; % frequencies to be tested
ARG.NPerm = 10000; % permutations for bootstrap analysis
ARG.Nsim = 1000; % number of simulation runs per parameters
ARG.Modeln = 5; % we assume 5 components in data simulation and analysis
ARG.Freq_U = 0.5; % frequency of U-term
ARG.Freq_Effect = 4;
ARG.Npad = 30;  % use 30 for DT of 60, or 45 for DT of 40

% -------------------------------------------------------------------
% simulate Nsub participants with effect at 2 Hz
% for each participajnt we simulate 750 trials with 
% gaussian distributed effect of rhythmicity

ARG.Nsub = 25;
ARG.Ntrial = 700; % assume we test each ear independently.
% mean, SD for each factor across the sample 
ARG.Factors(1,:) = [1,0.1]; % offset
ARG.Factors(2,:) = [0.1,0.1]; % slope
ARG.Factors(3,:) = [0.1,0.1]; % u/v
ARG.Factors(4,:) = [0.2,0.1]; % rythm
ARG.Factors(5,:) = [0.2,0.1]; % rythm

Nlev{1} = [2:8];
Nlev{2} = [2:2:6];

% -------------------------------------------------------
% set 1 - true effect, different noise levels

for n=1:length(Nlev{1})
 fprintf('set 1 %d \n',n);
 ARG.Noise = Nlev{1}(n);
  out={};
  parfor rep=1:ARG.Nsim
    out{rep} =  Binned_simulate(ARG);
  end
  Set1{n} = out;
end


% -------------------------------------------------------
% set 2 - absence of effect different noise levels
ARG.Factors(4,:) = [0.0,0.0];
ARG.Factors(5,:) = [0.0,0.0];
for n=1:length(Nlev{2})
  ARG.Noise = Nlev{2}(n);
  fprintf('set 2 %d \n',n);
  out={};
   parfor rep=1:ARG.Nsim
    out{rep} =  Binned_simulate(ARG);
  end
  Set2{n} = out;
end


save('New2Simulations_bins.mat','Set1','Set2','Nlev','ARG')


