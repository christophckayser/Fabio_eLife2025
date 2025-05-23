
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
Ind_TrialId = find(strcmp(ARG.VarList,'TotalTrialNr'));

% for dual task
Ind_task =  find(strcmp(ARG.VarList,'Timechange'));
Ind_stimchange = find(strcmp(ARG.VarList,'Fixchange'));
Ind_taskresp = find(strcmp(ARG.VarList,'ChangeResp'));


fieldname = {'dpSOA','critSOA','rtSOA'};
ARG.Do_shuffle = 0;
ARG.fast = 0;

% load  relevant Split information from previous run
sname = sprintf('%s/Timestats_eyesplit_Exp%d.mat',Prepropath,WHICH_EXP);
load(sname,'SPLITS');
Nsub = length(SPLITS.SplitTrials);

SubListOrig = SPLITS.isgood; % index in original sublist

% initialize variables
for SPLIT=1:2
  for task=1:2
    for par=1:3 % parameter
      for E=1:3 % each ear, or combined
        Vector{par,E,task} = zeros(Nsub,length(ARG.flist)+1);
        PvalMC{SPLIT}{par,E} = zeros(ARG.Nsim,length(ARG.flist));
      end
    end
  end
end
% in contrast to  Analyze_CK15_Exp4_eyesplits
% we here only use those n subjects retained in that function.
% speeds up things

for SPLIT=1:2 % which eye parameter
  % random participant sample
  for R=1:ARG.Nsim
    SubList = ceil(rand(1,Nsub)*Nsub); % random sample
    for S=1:Nsub
      Suse = SubList(S); % entry into SPLITS structure

      for task=1:2 % how to split the data
        if task==1
          usetrials = SPLITS.SplitTrials{Suse}{SPLIT}{1};
        else
          usetrials = SPLITS.SplitTrials{Suse}{SPLIT}{2};
        end
        % these are the actual trial IDs. Match to data
        TrialId = Data_all{SubListOrig(Suse)}(:,Ind_TrialId);
        [ia,ib,ic] = intersect(usetrials,TrialId);

        usetrials = ic;

        % ----------------------------------------------
        % VS for actual data
        Behav = compute_binned_behavior(Data_all{SubListOrig(Suse)}(usetrials,:),ARG);
        for par=1:3
          data = getfield(Behav,fieldname{par});
          for E=1:3 % both, left, right
            [~,vs] = local_fitmodelsGLM(data(E,:)',ARG);
            Vector{par,E,task}(S,:) = squeeze(sqrt(sum(vs([ARG.Modeln-1,ARG.Modeln],:).^2,1)));
          end % E
        end % par
      end % task

    end % S
    for par=1:3
      for E=1:3
        tmp = Vector{par,E,1} - Vector{par,E,2};
        df = size(tmp,1);
        pval = 2*tcdf(-abs(sqrt(df)*mean(tmp)./std(tmp)),df-1);
        PvalMC{SPLIT}{par,E}(R,:) = pval(2:end);
      end
    end

  end % R
end % SPLIT


sname = sprintf('%s/Timestats_eyesplitMC_Exp%d.mat',Prepropath,WHICH_EXP);
save(sname,'ARG','PvalMC')

return




