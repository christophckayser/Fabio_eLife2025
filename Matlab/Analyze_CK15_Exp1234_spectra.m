
clear;

% ---------------------------------------------------------------------
% Analysis of behavior time course  CK 15 Exp 1-4
% compute the spectrum and AR Surrogate data for each participant. 

% parameters
ARG.DT = 60; % time bin length to compute SOA binned metrics
ARG.BINS = [0:ARG.DT :1200]./1000;
ARG.NPerm = 10000; % number of permutations for surrogate data
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging acvros trials
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset
ARG.Npad = 30;  % use 30 for DT of 60, or 45 for DT of 40. Padding of frequency spectra

for WHICH_EXP = 1:4

  % --------------------------------------------------------------------------
  % load data
  defdirsCK15;
  sname = sprintf('%s/ProcessedData_Exp%d.mat',Prepropath,WHICH_EXP)
  x = load(sname);
  Data_all = x.Data_all;
  Nsub = length(Data_all);
  ARG.VarList = x.ARG.VarList;

  %----------------------------------------------------------------
  % initialize variables
  nf = (length(ARG.BINS)+ARG.Npad*2);
  if rem(nf,2)==1, nf = nf-1; end
  nf = nf/2-5;
  for par=1:3
    for E=1:3
      Spectra{par,E} = (zeros(Nsub,nf));
      SpectraShuf{par,E} = (zeros(Nsub,ARG.NPerm,nf));
    end
  end

  %----------------------------------------------------------------
  % compute spectra and surrogate
  fieldname = {'dpSOA','critSOA','rtSOA'};
  for S=1:Nsub
    Behav = compute_binned_behavior(Data_all{S},ARG);
    for par=1:3
      data = getfield(Behav,fieldname{par});
      for E=1:3
        [spec,fax,ARshuf] = local_computespectra(data(E,:)',ARG,ARG.NPerm);
        Spectra{par,E}(S,:) = spec;
        SpectraShuf{par,E}(S,:,[1:size(ARshuf,2)]) = ARshuf;
      end
    end
  end

  for par=1:3
    for E=1:3
      SpectraShuf{par,E}= single( SpectraShuf{par,E}); % save some memory
    end
  end

  sname = sprintf('%s/Spectra_Exp%d.mat',Prepropath,WHICH_EXP);
  save(sname,'ARG','SpectraShuf','Spectra','fax')
end


return;

