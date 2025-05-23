function Behav = local_compbehavsim(Data,ARG)

% compile single trial data into SOA-binned data
% potentially shuffle the data -> ARG.Do_shuffle
if ~isfield(ARG,'Do_shuffle')
  ARG.Do_shuffle= 0;
end

if ARG.Do_shuffle
  for s=1:length(Data)
    rp = randperm(size(Data{s},1));
    Data{s}(:,1) = Data{s}(rp,1)  ;
  end
end

ARG.nbins= length(ARG.BINS)-1;
Behav=[];
for s=1:length(Data)
  tmp = Data{s};
  for b=1:ARG.nbins
    % all trials with this SOA bin
    j = find( (tmp(:,2)>=ARG.BINS(b)).*(tmp(:,2)<ARG.BINS(b+1)));
    Behav(s,b) = mean(tmp(j,1));
  end
end
