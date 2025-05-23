
clear;

% ---------------------------------------------------------------------
% Analysis of behavior time course  CK 15 Exp 1-4
% USING BINNED DATA
% - fit the time course of SDT / RT against SOA using linear models
% - compare model fit and vector strength to test for evidence of rhythmicity
defdirsCK15;


% ----------------------------------------------------------------
% parameters
ARG.DT = 60; % SOA binning. each 60 ms, max frequency is 8.3 Hz
ARG.Modeln = 5; % number of model components in regression model.
ARG.flist = [1.2:0.1:4, 4.2:0.2:8]; % frequencies to be tested
ARG.NPerm = 10000; % number of permutations for surrogate data
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging acvros trials
ARG.Freq_U = 0.5; % frequency of U-term in regression model
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset

% SOA bins for analysis
ARG.BINS = [0:ARG.DT :1200]./1000;
tax = ARG.BINS(1:end-1)+0.03; % true time axis in experiment


% --------------------------------------------------------------------------

for WHICH_EXP = 1:4
  % load data

  sname = sprintf('%s/ProcessedData_Exp%d.mat',Prepropath,WHICH_EXP)
  x = load(sname);
  Data_all = x.Data_all;
  ARG.VarList = x.ARG.VarList;
  Nsub = length(Data_all);

 
  % --------------------------------------------------------------------------
  % 1) fit models with different temporal structures using GLMs
  % Fit_qual{parameter, analysis} (Sub,:,freq): [AIC];
  % Beta{param,analysis} (S,:,freq):  slopes
  % we do the analysis for both ears combined, and for each ear individually


  % initialize variables
  for par=1:3 % parameter
    for E=1:3 % each ear, or combined
      Fit_Qual{par,E} = zeros(Nsub,length(ARG.flist)+1);
      Vector{par,E} = zeros(Nsub,length(ARG.flist)+1);
      VectorShuf{par,E} = zeros(Nsub,ARG.NPerm,length(ARG.flist)+1);
    end
  end
  % AIC and Beta for actual data
  fieldname = {'dpSOA','critSOA','rtSOA'};
  ARG.Do_shuffle = 0;
  ARG.fast = 0;
  for S=1:Nsub
    Behav = compute_binned_behavior(Data_all{S},ARG);
    for par=1:3
      data = getfield(Behav,fieldname{par});
      for E=1:3 % both, left, right
        [Fit_Qual{par,E}(S,:,:),vs] = local_fitmodelsGLM(data(E,:)',ARG);
        Vector{par,E}(S,:) = squeeze(sqrt(sum(vs([ARG.Modeln-1,ARG.Modeln],:).^2,1)));
      end
    end
  end

  % Beta's for surrogate data
  ARG.Do_shuffle = 1;
  ARG.fast = 1;  % dont calculate LL
  for R=1:ARG.NPerm
    for S=1:Nsub
      Behavshuf = compute_binned_behavior(Data_all{S},ARG);
      for par=1:3
        data = getfield(Behavshuf,fieldname{par});
        for E=1:3
          [~,vs] = local_fitmodelsGLM(data(E,:)',ARG);
          VectorShuf{par,E}(S,R,:) = squeeze(sqrt(sum(vs([ARG.Modeln-1,ARG.Modeln],:).^2,1)));
        end
      end
    end
  end

  for par=1:3
    for E=1:3
      VectorShuf{par,E}= single( VectorShuf{par,E});
    end
  end
  sname = sprintf('%s/Binned_Exp%d.mat',Prepropath,WHICH_EXP);
  save(sname,'Fit_Qual','ARG','Vector','VectorShuf')
end


return;

