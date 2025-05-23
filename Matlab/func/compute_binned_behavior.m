
function Behav = compute_binned_behavior(data,ARG)


% function Behav = compute_binned_behavior(data,ARG)
% 
% summarize behavioral data for binned SOAs
% Compute SDT and RT per ear and for both ears
% 
% ARG.BINS
% ARG.rt_transform
% ARG.Do_shuffle : 1 for shuffling SOA ACROSS TRIALS. Used only for stats

nbins = length(ARG.BINS)-1;
% --------------------------------------------------------------------------
% data format for Exp 1/2
% Result(trial,:) = [trial,condMAT(trial,1:5),RT,resp,resp_code,Tone_SNR,trial,subid,RT_Target];
% CondMat  [Freq, Ear, SOA, SOABin, seed left, seed right]
%  -> Freq (1 Low, 2 High), Ear (1 Left, 2 Right), SOA

% data format for Exp 3
% Result(trial,:) = [trial,condMAT(trial,1:5),RT,resp,resp_code,Freq_diff,,trial,subid,RT_Target];
% CondMat  [Freq, Ear, SOA, SOABin, seed left, seed right]
%  -> Freq (1 st is low, 2 is low), Ear (1 Left, 2 Right), SOA

% exp 4

% Result(trial,:) = [trial,condMAT(trial,1:5),RT,resp,resp_code,Tone_SNR,fix_change,tchange,respchange];
% Timing(trial,:) = [t1,t2,t3,t4,StartTime,t5];

% set indices of relevant variables
Ind_Freq =  find(strcmp(ARG.VarList,'Freq'));
Ind_Ear =  find(strcmp(ARG.VarList,'Ear'));
Ind_SOA =  find(strcmp(ARG.VarList,'SOA'));
Ind_SOABin =  find(strcmp(ARG.VarList,'SOABin'));
Ind_Resp = find(strcmp(ARG.VarList,'Resp'));
Ind_Thr = find(strcmp(ARG.VarList,'SNR'));

Ind_RT =  find(strcmp(ARG.VarList,ARG.WHICHRT));

% --------------------------------------------------------------------------
% Shuffle SOA if required
% --------------------------------------------------------------------------
if ~isfield(ARG,'Do_shuffle')
  ARG.Do_shuffle = 0;
end
if ARG.Do_shuffle
  rp = randperm(size(data,1));
  data(:,[Ind_SOA,Ind_SOABin]) = data(rp,[Ind_SOA,Ind_SOABin])  ;
end

% --------------------------------------------------------------------------
% overall performance independent of SOA
% --------------------------------------------------------------------------
RTS = data(:,Ind_RT);
RTS = ARG.rt_transform(RTS);
Perf = data(:,Ind_Resp); % correct, wrong response

for E=1:2 % each EAR
  for f=1:2
    % 1) pc per frequency and Ear
    t = find( (data(:,Ind_Freq)==f).*(data(:,Ind_Ear)==E));
    Behav.pc(f,E) = mean(Perf(t));
    % rt for correct responses
    tmp = RTS(t);
    Behav.rt(f,E) = mean(tmp);
  end

  % 2) signal detection theory per Ear
  % we compute hit rate based on f==1
  t1 = find( (data(:,Ind_Freq)==1).*(data(:,Ind_Ear)==E));
  HR = mean(Perf(t1)==1); % hit rate for f=1
  t2 = find( (data(:,Ind_Freq)==2).*(data(:,Ind_Ear)==E));
  FR =  mean(Perf(t2)==0); % false alarm rate for f=1.
  ntrl = (length(t1)+length(t2))/2;
  [d,beta,c,A,B] = ck_stat_signalTheory(HR,FR,ntrl);
  Behav.dp(E) = d;
  Behav.crit(E) = c;
end

% --------------------------------------------------------------------------
% performance as a function of SOA - both ears combined
% --------------------------------------------------------------------------

for b=1:nbins
  % all trials with this SOA bin
  j = find( (data(:,Ind_SOA)>=ARG.BINS(b)).*(data(:,Ind_SOA)<ARG.BINS(b+1)));
  tmp = data(j,:);

  % signal detection theory
  t1 = find( (tmp(:,Ind_Freq)==1));
  HR = nanmean(Perf(j(t1))==1); % hit rate for f=1
  t2 = find( (tmp(:,Ind_Freq)==2));
  FR =  nanmean(Perf(j(t2))==0); % false alarm rate for f=1
  ntrl = (length(t1)+length(t2))/2;
  [d,beta,c,A,B] = ck_stat_signalTheory(HR,FR,ntrl);
  Behav.dpSOA(1,b) = d;
  Behav.critSOA(1,b) = c;

  % reaction time and performance
  tmprt = RTS(j);
  Behav.rtSOA(1,b) = mean(tmprt);
  Behav.ntrial(1,b) = length(j);
  Behav.PC(1,b) = nanmean(Perf(j));
end

% --------------------------------------------------------------------------
% performance as a function of SOA - each ear
% --------------------------------------------------------------------------
bins = unique(data(:,Ind_SOABin));
for E=1:2
  for b=1:nbins
    % all trials with this SOA bin
    j = find( (data(:,Ind_SOA)>=ARG.BINS(b)).*(data(:,Ind_SOA)<ARG.BINS(b+1)).*(data(:,Ind_Ear)==E));
    tmp = data(j,:);

    % signal detection theory
    t1 = find( (tmp(:,Ind_Freq)==1));
    HR = nanmean(Perf(j(t1))==1); % hit rate for f=1
    t2 = find( (tmp(:,Ind_Freq)==2));
    FR =  nanmean(Perf(j(t2))==0); % false alarm rate for f=1
    ntrl = (length(t1)+length(t2))/2;
    [d,beta,c,A,B] = ck_stat_signalTheory(HR,FR,ntrl);
    Behav.dpSOA(E+1,b) = d;
    Behav.critSOA(E+1,b) = c;

    % reaction time and performance
    tmprt = RTS(j);
    Behav.rtSOA(E+1,b) = mean(tmprt);
    Behav.ntrial(E+1,b) = length(j);
    Behav.PC(E+1,b) = nanmean(Perf(j));

  end
end
% to store the effective thresholds, we bin them into 10 time bins over
% the course of all included trials
tmp = data(:,Ind_Thr);  n = floor(size(tmp,1)/10);
tmp = reshape(tmp([1:n*10]),[n,10]);

BehavSoa.Thres = mean(tmp,1);

