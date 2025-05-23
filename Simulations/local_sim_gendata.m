function Data = local_sim_gendata(ARG)

% simulate singel trial RTs for ARG.Nsub participants based on a linear model
% the parameters are coded as ARG.Factors(1:4, [mean, SD]) 1:5 as in the
% linear model
% add measurement noise ARG.Noise

% use actual condition matrix
Cond = load('CondSim.mat');
for s=1:ARG.Nsub
  % draw beta's
  beta=[];
  for b=1:5
    tmp = randn(100,1)*ARG.Factors(b,2)+ARG.Factors(b,1);
    beta(b) = tmp(1);
  end
  % we simulate data based on this model, given the trial specific SOAbin
  soa = Cond.Conditions{s}([1:ARG.Ntrial],3);
  X = zeros(ARG.Ntrial,5);
  X(:,1) = 1;
  X(:,2) = soa; % linear
  X(:,3) = cos(2*pi*soa*ARG.Freq_U);
  X(:,4) = cos(2*pi*soa*ARG.Freq_Effect);
  X(:,5) = sin(2*pi*soa*ARG.Freq_Effect);

  Data{s}(:,:) = [ X*beta',soa, soa];

  % add measurement noise
  Data{s}(:,1) = Data{s}(:,1) +randn(size(Data{s}(:,2) )) * ARG.Noise;
end