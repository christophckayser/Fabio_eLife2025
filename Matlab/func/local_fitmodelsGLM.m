function [AIC, beta] = local_fitmodelsGLM(y,ARG)

% function [AIC, beta] = local_fitmodelsGLM(y,ARG)
%
% fits linear models to the BINNED behavioral data and returns slopes and AIC
%
% input: 
%    data  expected  [RT, SOA, time_axis(sec)])
%    ARG.flist     frequencies to be tested
%    ARG.Modeln    Number of components for model. 4) or 5) to also include u/v shaped 
%    ARG.fast    1) to only compute beta, 0) to also compute AIC
% 
% out: 
%    AIC(freq) and  Beta(:,freq) freq=1 constant model


n = length(y);
nf = length(ARG.flist);
beta = zeros(ARG.Modeln,nf+1);
AIC = zeros(nf+1,1);
tax = ARG.BINS(1:end-1)+0.03;

% -----------------------------------------------------------------------------
% fit 'trivial models'. These explain the data based on a 
% - constant + linear slope + a U/V shaped term adapted to the stimulus duration

if ARG.Modeln==4
  X(:,1) = ones(n,1);
  X(:,2) = tax; % linear time
elseif ARG.Modeln==5
  X(:,1) = ones(n,1);
  X(:,2) = tax; % linear
  X(:,3) = cos(2*pi*tax*ARG.Freq_U ); % u/v
end
k = ARG.Modeln-2;
[beta([1:ARG.Modeln-2],1),ll]  = fast_regress(y,X,ARG.fast);
AIC(1) =  -2*ll(1) + 2*k;% + 2*k*(k+1)/(n-k-1);

% models with rhytymic component
k = ARG.Modeln;
for f=1:nf
  X(:,ARG.Modeln-1) = cos(2*pi*tax*ARG.flist(f));
  X(:,ARG.Modeln)   = sin(2*pi*tax*ARG.flist(f));
  [beta(:,f+1),ll]= fast_regress(y,X,ARG.fast);
  AIC(f+1) =  -2*ll(1) + 2*k;% + 2*k*(k+1)/(n-k-1);
end


return;










