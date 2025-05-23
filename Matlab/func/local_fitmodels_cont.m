function [AIC, beta] = local_fitmodels_cont(data,ARG)

% function [AIC, beta] = local_fitmodels_cont(data,ARG)
%
% fits linear models to the continous data and returns slopes and AIC
%
% input: 
%    data   (time x 2) with dim 1 data dim 2 time in sec
%    ARG.flist     frequencies to be tested
%    ARG.Modeln    Number of components for model. 4) or 5) to also include u/v shaped 
%    ARG.fast    1) to only compute beta, 0) to also compute AIC
% 
% out: 
%    AIC(freq) and  Beta(:,freq) freq=1 constant model


nf = length(ARG.flist);
beta = zeros(ARG.Modeln,nf+1);
AIC = zeros(nf+1,1);
n = size(data,1);

if ARG.Modeln==4
  X = ones(n,1);
  X(:,2) = data(:,2); % linear time
elseif ARG.Modeln==5
  X = ones(n,1);
  X(:,2) = data(:,2); % linear time
  X(:,3) = cos(2*pi*data(:,2)*ARG.Freq_U);
end
k = ARG.Modeln-2;
[beta([1:ARG.Modeln-2],1),ll] = fast_regress(data(:,1),X,ARG.fast);
AIC(1) =  -2*ll(1) + 2*k;% + 2*k*(k+1)/(n-k-1);

% models with rhytymic component
k = ARG.Modeln;
for f=1:nf
  X(:,ARG.Modeln-1) = cos(2*pi*data(:,2)*ARG.flist(f));
  X(:,ARG.Modeln)   = sin(2*pi*data(:,2)*ARG.flist(f));
  [beta(:,f+1),ll] = fast_regress(data(:,1),X,ARG.fast);
  AIC(f+1) = -2*ll(1)  + 2*k ;%+ 2*k*(k+1)/(n-k-1);
end


return;












