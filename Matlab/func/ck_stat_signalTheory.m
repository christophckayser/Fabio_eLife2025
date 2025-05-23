function [d,beta,c,A,B] = ck_stat_signalTheory(HitR,FaR,Ntrials)
% Signal Detection Theory
%
% function [d,beta,c,A,B] = ck_stat_signalTheory(HitR,FaR)
%
% signal detection theory analysis. HitR = Hitrate (in proportion)
%
% hit rate H: proportion of STIM trials to which subject responded YES =  P("yes" | YES) 
% false alarm rate FaR: proportion of NO-STIM trials to which subject responded YES  = P("yes" | NO)
%
% output: d' (sensitivity), beta (bias), c (criterion),
%         A (non-parametric sensitvity), B'' (non-parametric bias)
%
%
% Example:               stim (Y/N)
%                    80        60  
%          resp(Y/N) 20        40
%
% -> HR=  80/100=0.8    FR = 60/100=0.6


% compute data augmentation term 
if nargin > 2
    AugT = 1/(Ntrials * 2);
else
    AugT = eps;
end

% d'
if HitR==1
  HitR = HitR-AugT;
elseif HitR==0
  HitR = HitR+AugT;
end
if FaR==1
  FaR = FaR-AugT;
elseif FaR==0
  FaR = FaR+AugT;
end
d = norminv(HitR) - norminv(FaR);

% A'
A = 0.5 + sign(HitR-FaR) * ( (HitR-FaR)^2 + abs(HitR-FaR))/(4*max(HitR,FaR)-4*HitR*FaR);

% beta
beta = exp( ((norminv(FaR)^2) - (norminv(HitR)^2))/2 );

% criterion 
c = - (norminv(HitR) + norminv(FaR)) / 2;

% B''
B = sign(HitR-FaR) * ( (HitR*(1-HitR) - FaR*(1-FaR)) / (HitR*(1-HitR) + FaR*(1-FaR)));



