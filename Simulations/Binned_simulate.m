
function out =  Binned_simulate(ARG)

% simulates data and implements diffrent analyses
%
% output:
%  GLM results:
%    out{1} AIC
%    out{2} PvalBeta;
%  Spectra results:
%    out{3} = PvalSpecAR;
%    out{4} = fax_spec;


% -------------------------------------
% generate data for n participants
Data = local_sim_gendata(ARG);



% ------------------------------------------------------------------
% fit  models

Beta = zeros(ARG.Nsub,ARG.Modeln,length(ARG.flist)+1);
BetaShuf =  zeros(ARG.Nsub,ARG.Modeln,length(ARG.flist)+1,ARG.NPerm);
Fit_Qual = zeros(ARG.Nsub,length(ARG.flist)+1);

% convert to time-binned data
ARG.Do_shuffle=0;
Behav = local_compbehavsim(Data,ARG);
ARG.fast = 0;
for s=1:ARG.Nsub
  [Fit_Qual(s,:),Beta(s,:,:)] = local_fitmodelsGLM(Behav(s,:)',ARG);
end
out{1}  = sq(sum(Fit_Qual)); % group-level fit


% Beta for surrogate data
ARG.Do_shuffle = 1;
ARG.fast = 1;
for rep=1:ARG.NPerm
  BehavS = local_compbehavsim(Data,ARG);
  for s=1:ARG.Nsub
    [~,BetaShuf(s,:,:,rep)] = local_fitmodelsGLM(BehavS(s,:)',ARG);
  end
end

% p-value for vector strength
Vector = squeeze(sqrt(sum(Beta(:,[ARG.Modeln-1,ARG.Modeln],[2:end]).^2,2)));
VectorShuf = squeeze(sqrt(sum(BetaShuf(:,[ARG.Modeln-1,ARG.Modeln],[2:end],:).^2,2)));
Vector = squeeze(mean(Vector,1));
VectorShuf = squeeze(mean(VectorShuf,1));

PvalBeta=[];
for f=1:length(ARG.flist)
  PvalBeta(f) = mean(VectorShuf(f,:)>Vector(f));
end
out{2} = PvalBeta;


% -------------------------------------------------------
%  AR spectrum based on time-binned data
ARG.Do_shuffle=0;
% actual data
data = local_compbehavsim(Data,ARG);
[spec,fax_spec] = local_computespectra(data(1,:)',ARG); % dummy run
nf = length(fax_spec);

Spectra = zeros(ARG.Nsub,nf);
SpectraAR= zeros(ARG.Nsub,ARG.NPerm,nf);
for S=1:ARG.Nsub
  [Spectra(S,:) ,~,SpectraAR(S,:,:)] = local_computespectra(data(S,:)',ARG,ARG.NPerm);
end

% p - value spectra
Spectrue = squeeze(mean(Spectra,1));
SpecShuf  = squeeze(mean(SpectraAR,1));
PvalSpecAR=[];
for f=1:length(fax_spec)
  PvalSpecAR(f) = mean(SpecShuf(:,f)>Spectrue(f));
end


out{3} = PvalSpecAR;
out{4} = fax_spec;
