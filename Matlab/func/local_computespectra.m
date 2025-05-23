function [spec,fax,PowerAR] = local_computespectra(y,ARG,doAR)

% function [spec,fax,PowerAR] = local_computespectra(y,ARG,doAR)
%
% compute frequency spectrum for behavioral data. 
%
% input: y behavioral data per SOA bin
% doAR > 0 - implement ARbased estiamte of spectral estimate using doAR simulations
% in line with other analyses, we retain only f>1Hz

if nargin < 3
  doAR = 0;
end
y = y(:)';
npad = ARG.Npad;


% standardize overall power
tmp = detrend(y);
tmp = tmp-mean(tmp(:));
% frequency spectrum after zero-padding
[spec,fax] = ck_powerspec([zeros(1,npad),tmp,zeros(1,npad)],1000./ARG.DT);
% DC and <1 Hz removed 
keep = find(fax>1.0);
spec = spec(keep);
fax = fax(keep);


PowerAR=[];
% ----------------------------------------
% surrogate data. requires AR Toolbox
if doAR
  PowerAR = zeros(doAR,length(fax));
  % fit using the ARfit toolbox
  [w,A,C]=arfit(tmp',1,1) ;
  % sim model
  for R=1:doAR
    % genreate data using this model
    [v] = arsim(w,A,C,length(tmp),2000);
    % frequency spectrum after zero-padding to fixed common resolution across paradigms
    [p,~] = ck_powerspec([zeros(1,npad),v',zeros(1,npad)],1000./ARG.DT);
    PowerAR(R,:) = p(keep);
  end
end

return