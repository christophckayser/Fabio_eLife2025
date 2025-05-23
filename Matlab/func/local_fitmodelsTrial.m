function [AIC, beta] = local_fitmodelsTrial(data,ARG)

% function [AIC, beta] = local_fitmodelsGLM(data,ARG)
%
% fits linear models to the SINGEL TRIAL data and returns slopes and AIC
% ARG.Modeltype==2 for RT, 1 for Accuracy
% input: 
%    data  expected  [RT, Accuracy, SOA])
%    ARG.flist     frequencies to be tested
%    ARG.Modeln    Number of components for model. 4) or 5) to also include u/v shaped 
%    ARG.Modeltype  1) linar model for rt 2) logistic model for accuracy
%    ARG.fast    1) to only compute beta, 0) to also compute AIC
% 
% out: 
%    AIC(freq) and  Beta(:,freq) freq=1 constant model


n = size(data,1);
nf = length(ARG.flist);
beta = zeros(ARG.Modeln,nf+1);
AIC = zeros(nf+1,1);


if ARG.Modeltype==2
  % -------------------------------------------------------------
  % linear models for RT
  % -------------------------------------------------------------
  if ARG.Modeln==4
    X = ones(n,1);
    X(:,2) = data(:,3); % linear time
 elseif ARG.Modeln==5
    X = ones(n,1);
    X(:,2) = data(:,3); % linear time
    X(:,3) = cos(2*pi*data(:,3)*ARG.Freq_U);
  end
  k = ARG.Modeln-2;
  [beta([1:ARG.Modeln-2],1),ll] = fast_regress(data(:,1),X,ARG.fast);
  AIC(1) =  -2*ll(1) + 2*k;% + 2*k*(k+1)/(n-k-1);

  % models with rhytymic component
  k = ARG.Modeln;

  for f=1:nf
    X(:,ARG.Modeln-1) = cos(2*pi*data(:,3)*ARG.flist(f));
    X(:,ARG.Modeln)   = sin(2*pi*data(:,3)*ARG.flist(f));
    [beta(:,f+1),ll] = fast_regress(data(:,1),X,ARG.fast);
    AIC(f+1) = -2*ll(1)  + 2*k ;%+ 2*k*(k+1)/(n-k-1);
  end

elseif ARG.Modeltype==1

  % ------------------------------------------
  % linear models for accuracy
  if ARG.Modeln==4
    X(:,1) = data(:,3); % linear time
  elseif ARG.Modeln==5
    X(:,1) = data(:,3); % linear time
    X(:,2) = cos(2*pi*data(:,3)*ARG.Freq_U);
  end
  k = ARG.Modeln-2;
  [b,~,ll] = ck_decoding_logist(X,data(:,2),[],0,0);
  beta([1:length(b)],1) = b([end,1:end-1]); % because constant is put last
  AIC(1) =  -2*ll(1) + 2*k ;% + 2*k*(k+1)/(n-k-1);

  % models with rhytymic component
  k = ARG.Modeln;
  for f=1:nf
    X(:,ARG.Modeln-2) = cos(2*pi*data(:,3)*ARG.flist(f));
    X(:,ARG.Modeln-1) = sin(2*pi*data(:,3)*ARG.flist(f));
    [b,~,ll] =  ck_decoding_logist(X,data(:,2),[],0,0);
    beta([1:length(b)],f+1) = b([end,1:end-1]); % because constant is put last
    AIC(f+1) = -2*ll(1)  + 2*k  ;%+ 2*k*(k+1)/(n-k-1);
  end

end



return;












