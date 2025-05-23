
clear;

% ---------------------------------------------------------------------
% Monte Carlo simulations of participant samples


% parameters
ARG.DT = 60; % time bin length.
ARG.BINS = [0:ARG.DT :1200]./1000;
ARG.NPerm = 10000; % number of permutations for surrogate data
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging acvros trials
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset
ARG.Npad = 30;  % use 30 for DT of 60, or 45 for DT of 40

ARG.Nsim = 5000; % number of draws of participant samples

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
      Spectra{par,E} = zeros(Nsub,nf);
      SpectraMC{par,E} = zeros(ARG.Nsim,nf);
    end
  end

  %----------------------------------------------------------------
  % compute spectra 
  fieldname = {'dpSOA','critSOA','rtSOA'};
  for R=1:ARG.Nsim
    SubList = ceil(rand(1,Nsub)*Nsub);
    for S=1:Nsub
      Suse = SubList(S);
      Behav = compute_binned_behavior(Data_all{Suse},ARG);
      for par=1:3
        data = getfield(Behav,fieldname{par});
        for E=1:3
          [spec,fax,~] = local_computespectra(data(E,:)',ARG);
          Spectra{par,E}(S,:) = spec;
        end
      end
    end
    for par=1:3
      for E=1:3
        SpectraMC{par,E}(R,:) = mean(Spectra{par,E});
      end
    end
  end

  sname = sprintf('%s/SpectraMC_Exp%d.mat',Prepropath,WHICH_EXP);
  save(sname,'ARG','SpectraMC','fax')
end


return;

