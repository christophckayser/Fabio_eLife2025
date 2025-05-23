
clear;
% ---------------------------------------------------------------------
% Analysis of behavior time course for CK 15 Exp 4
% comparison of both experimental conditions
% rathern than comparing to a surrogate, we compare the 
% effect strength between conditions

WHICH_EXP = 4;
ARG.DT = 60; % each 60 ms, max frequency is 8.3 Hz

ARG.Modeln = 5; % number of model components
ARG.flist = [1.2:0.1:4, 4.2:0.2:8]; % frequencies to be tested
ARG.NPerm = 10000; % number of permutations for surrogate data
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging acvros trials
ARG.Freq_U = 0.5; % frequency of U-term
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset

ARG.BINS = [0:ARG.DT :1200]./1000;
tax = ARG.BINS(1:end-1)+0.03; % true time axis in experiment
ARG.Nsim = 5000;

% --------------------------------------------------------------------------
% load data
defdirsCK15;
sname = sprintf('%s/ProcessedData_Exp%d.mat',Prepropath,WHICH_EXP);
x = load(sname);
Data_all = x.Data_all;
ARG.VarList = x.ARG.VarList;

Ind_Freq =  find(strcmp(ARG.VarList,'Freq'));
Ind_Ear =  find(strcmp(ARG.VarList,'Ear'));
Ind_SOA =  find(strcmp(ARG.VarList,'SOA'));
Ind_SOABin =  find(strcmp(ARG.VarList,'SOABin'));
Ind_Resp = find(strcmp(ARG.VarList,'Resp'));
Ind_Thr = find(strcmp(ARG.VarList,'SNR'));
Ind_RT =  find(strcmp(ARG.VarList,ARG.WHICHRT));

% for dual task
Ind_task =  find(strcmp(ARG.VarList,'Timechange'));
Ind_stimchange = find(strcmp(ARG.VarList,'Fixchange'));
Ind_taskresp = find(strcmp(ARG.VarList,'ChangeResp'));

Nsub = length(Data_all);

fieldname = {'dpSOA','critSOA','rtSOA'};
ARG.Do_shuffle = 0;
ARG.fast = 1; %1

% initialize variables
for par=1:2 % parameter
  for E=1:3 % each ear, or combined
    for task=1:2
      Vector{par,E,task} = zeros(Nsub,length(ARG.flist)+1);
    end
    PvalMC{par,E} = zeros(ARG.Nsim,length(ARG.flist));
  end
end



% random participant sample
for R=1:ARG.Nsim
  SubList = ceil(rand(1,Nsub)*Nsub);
  for S=1:Nsub
    Suse = SubList(S);
    for task=1:2
      if task==1 % with task
        usetrials = find(~isnan(Data_all{Suse}(:,Ind_task)));
      else % no task
        usetrials = find(isnan(Data_all{Suse}(:,Ind_task)));
      end
      % ----------------------------------------------
      % VS for actual data
      Behav = compute_binned_behavior(Data_all{Suse}(usetrials,:),ARG);
      for par=1:3
        data = getfield(Behav,fieldname{par});
        for E=1:3 % both, left, right
          [~,vs] = local_fitmodelsGLM(data(E,:)',ARG);
          Vector{par,E,task}(S,:) = squeeze(sqrt(sum(vs([ARG.Modeln-1,ARG.Modeln],:).^2,1)));
        end
      end % par
    end % task
  end % S


  % we need the t-value between tasks
  for par=1:3
    for E=1:3
      tmp = Vector{par,E,1} - Vector{par,E,2};
      df = size(tmp,1);
      pval = 2*tcdf(-abs(sqrt(df)*mean(tmp)./std(tmp)),df-1);
      PvalMC{par,E}(R,:) = pval(2:end);
    end
  end
end % R




sname = sprintf('%s/Timestats_ConditionsMC_Exp%d.mat',Prepropath,WHICH_EXP);
save(sname,'ARG','PvalMC')



return

