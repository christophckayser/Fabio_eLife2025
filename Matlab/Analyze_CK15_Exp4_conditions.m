
clear;
% ---------------------------------------------------------------------
% Analysis of behavior time course for CK 15 Exp 4
% comparison of both experimental conditions
% rathern than comparing to a surrogate, we compare the 
% effect strength between conditions
% uses analysis of binned SOA data

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
ARG.fast = 0;

% initialize variables
for task=1:2
  for par=1:2 % parameter
    for E=1:3 % each ear, or combined
      Vector{task}{par,E} = zeros(Nsub,length(ARG.flist)+1);
    end
  end
end


for S=1:Nsub
  for task=1:2
    if task==1 % with task
      usetrials = find(~isnan(Data_all{S}(:,Ind_task)));
    else % no task
      usetrials = find(isnan(Data_all{S}(:,Ind_task)));
    end

    % ----------------------------------------------
    % VS for actual data
    Behav = compute_binned_behavior(Data_all{S}(usetrials,:),ARG);
    for par=1:3
      data = getfield(Behav,fieldname{par});
      for E=1:3 % both, left, right
        [~,vs] = local_fitmodelsGLM(data(E,:)',ARG);
        Vector{task}{par,E}(S,:) = squeeze(sqrt(sum(vs([ARG.Modeln-1,ARG.Modeln],:).^2,1)));
        Ntrials(S,task) = length(usetrials);

      end % E
    end

    % -------------------------------------------------
    % compare auditory task performance
    Performance_aud.dp(S,task) = mean(Behav.dp);
    Performance_aud.crit(S,task) = mean(Behav.crit);
    Performance_aud.rt(S,task) = mean(Behav.rt(:));
  
  end % task



  % -------------------------------------------------
  % compute behavioral performance for change task

  trls = ~isnan(Data_all{S}(:,Ind_task));
  stim = Data_all{S}(trls,Ind_stimchange);
  % for trials with stim ~=0 we expect key 38, for those with 0 key 40
  resp = Data_all{S}(trls,Ind_taskresp);
  HR = sum(resp(stim~=0)==38)/sum(stim~=0);
  FR = sum(resp(stim==0)==38)/sum(stim==0);
  [d,~,c] = ck_stat_signalTheory(HR,FR,length(resp));

  Performance_task(S,:) = [d,c];



end  % S




sname = sprintf('%s/Timestats_Conditions_Exp%d.mat',Prepropath,WHICH_EXP);
save(sname,'ARG','Vector','Performance_task')




% reporting of task performance
[h,p,~,stat] = ttest(Performance_aud.dp(:,1),Performance_aud.dp(:,2));
fprintf('D-prime %1.2f+%1.2f vs. %1.2f+%1.2f, t=%1.2f, p=%1.4f \n', ... 
  mean( Performance_aud.dp(:,1)), sem( Performance_aud.dp(:,1)), mean(Performance_aud.dp(:,2)), sem(Performance_aud.dp(:,2)), stat.tstat,p);

[h,p,~,stat] = ttest(Performance_aud.crit(:,1),Performance_aud.crit(:,2));
fprintf('Bias %1.2f+%1.2f vs. %1.2f+%1.2f, t=%1.2f, p=%1.4f \n', ... 
  mean( Performance_aud.crit(:,1)), sem( Performance_aud.crit(:,1)), mean(Performance_aud.crit(:,2)), sem(Performance_aud.crit(:,2)), stat.tstat,p);

[h,p,~,stat] = ttest(Performance_aud.rt(:,1),Performance_aud.rt(:,2));
fprintf('Rt %1.3f+%1.3f vs. %1.3f+%1.3f, t=%1.2f, p=%1.4f \n', ... 
  mean( Performance_aud.rt(:,1)), sem( Performance_aud.rt(:,1)), mean(Performance_aud.rt(:,2)), sem(Performance_aud.rt(:,2)), stat.tstat,p);


return






