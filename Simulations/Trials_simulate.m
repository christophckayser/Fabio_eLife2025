
function out =  Trials_simulate(ARG)

% simulates data and implements trial-wise analysis
%
% output:
%  GLM results:
%    out{1} AIC
%    out{2} PvalBeta;


% -------------------------------------------------------------
% generate data for n participants
Data = local_sim_gendata(ARG);


% --------------------------------------------------------------------------
% compute model fit and beta for the actual data
Fit_Qual = zeros(ARG.Nsub ,length(ARG.flist)+1);
Beta = zeros(ARG.Nsub,ARG.Modeln,length(ARG.flist)+1);
BetaShuf =  zeros(ARG.Nsub,ARG.Modeln,length(ARG.flist)+1,ARG.NPerm);
for S=1:ARG.Nsub
  tmp = Data{S};
  ARG.fast = 0;
  [Fit_Qual(S,:),Beta(S,:,:)] = local_fitmodelsTrial(tmp,ARG);

  % Beta for surrogate data
  ARG.fast = 1;
  for rep=1:ARG.NPerm
    R = randperm(size(tmp,1));
    tmp(:,[1,2]) = tmp(R,[1,2]);
    [~,BetaShuf(S,:,:,rep)] = local_fitmodelsTrial(tmp,ARG);
  end
end
out{1}  = sq(sum(Fit_Qual)); % group-level fit

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

