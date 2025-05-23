
clear;
% ---------------------------------------------------------------------
% Analysis of behavior time course for CK 15 Exp 4
%
% separate fit to trials selected based on eye movement properties
% we restrict to participants with sufficient eye movement trials 
% SPLIT 1 Pupil dilation
% SPLIT 2 position SD


WHICH_EXP = 4;

ARG.Modeln = 5; % number of model components
ARG.flist = [1.2:0.1:4, 4.2:0.2:8]; % frequencies to be tested
ARG.NPerm = 8000; % number of permutations for surrogate data
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging acvros trials
ARG.Freq_U = 0.5; % frequency of U-term
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset
ARG.DT = 60; % each 60 ms, max frequency is 8.3 Hz
ARG.BINS = [0:ARG.DT :1200]./1000;
tax = ARG.BINS(1:end-1)+0.03; % true time axis in experiment

% SubList = {'A01','A02','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A32','A33','A34','A35','A36','A37'};

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
Ind_TrialId = find(strcmp(ARG.VarList,'TotalTrialNr'));

% for dual task
Ind_task =  find(strcmp(ARG.VarList,'Timechange'));
Ind_stimchange = find(strcmp(ARG.VarList,'Fixchange'));
Ind_taskresp = find(strcmp(ARG.VarList,'ChangeResp'));


fieldname = {'dpSOA','critSOA','rtSOA'};
ARG.Do_shuffle = 0;
ARG.fast = 0;

% save split trials for separate analysis
sname = sprintf('%s/Split_trials_Exp4.mat',Prepropath);
SPLITS = load(sname);

Nsub = length(SPLITS.SubListUse);
SubUse = SPLITS.isgood; % index into the original sublist


% initialize variables
for SPLIT=1:2
  for task=1:2
    for par=1:2 % parameter
      for E=1:3 % each ear, or combined
        Vector{task,SPLIT}{par,E} = zeros(20,length(ARG.flist)+1);
      end
    end
  end
end

Nsplit = [1,2];
for SPLIT=1:2 % which eye factor to use 
  for S=1:Nsub
    for task=1:2 % how to split the data
      if task==1
        usetrials = SPLITS.SplitTrials{S}{SPLIT}{1};
        % these are the actual trial IDs. Match to data
      else
        usetrials = SPLITS.SplitTrials{S}{SPLIT}{2};
      end
      TrialId = Data_all{SubUse(S)}(:,Ind_TrialId); 
      [ia,ib,ic] = intersect(usetrials,TrialId);

      usetrials = ic;
      Ntrials(S,SPLIT,task) = length(usetrials);

      % ----------------------------------------------
      % VS for actual data
      Behav = compute_binned_behavior(Data_all{SubUse(S)}(usetrials,:),ARG);
      for par=1:3
        data = getfield(Behav,fieldname{par});
        for E=1:3 % both, left, right
          [~,vs] = local_fitmodelsGLM(data(E,:)',ARG);
          Vector{task,SPLIT}{par,E}(S,:) = squeeze(sqrt(sum(vs([ARG.Modeln-1,ARG.Modeln],:).^2,1)));        
        end % E
      end
    end % task
  end  % S
end % split


% 
% find participants with sufficient many trials per Eear and split
crit = min(min(Ntrials,[],3),[],2);

Ncrit = 400;
isbad = find( crit<Ncrit);
isgood =  find( crit>Ncrit);


for SPLIT=1:2
  for task=1:2
    for par=1:3
      for E=1:3
        Vector{task,SPLIT}{par,E}(isbad,:) = [];
      end
    end
  end
end

% reduce SPLIT structure for use in MC simulations
SPLITS.isgood = SPLITS.isgood(isgood)
SPLITS.SubListUse = {SPLITS.SubListUse{isgood}};
SPLITS.SplitTrials = {SPLITS.SplitTrials{isgood}};

sname = sprintf('%s/Timestats_eyesplit_Exp%d.mat',Prepropath,WHICH_EXP);
save(sname,'ARG','Vector','SPLITS')

return








% OLD CODE FOR SINGEL TRIAL
% 
% 
%       tmp = Data_all{S}(usetrials,[Ind_RT,Ind_Resp,Ind_SOA,Ind_Thr,Ind_Ear]);
%       tmp(:,1) = ARG.rt_transform(tmp(:,1));
%       
%       for E=1:3
%         if E<3
%           j = find( (tmp(:,5)==E));
%         else
%           j = [1:size(tmp,1)];
%         end
%         for par=1:2
%           ARG.Modeltype = par; ARG.fast = 0;
%           [Fit_Qual{task,SPLIT}{par,E}(S,:),vs] = local_fitmodelsTrial(tmp(j,[1:3]),ARG);
%           Vector{task,SPLIT}{par,E}(S,:) = squeeze(sqrt(sum(vs([ARG.Modeln-1,ARG.Modeln],:).^2,1)));
%         end
%         Ntrials(S,E,SPLIT,task) = length(j);
