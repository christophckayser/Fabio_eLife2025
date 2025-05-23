
clear;

% CK 15 Exp 4
% Properties of eye movements in trial prior to target
% collects eye movement data and saves file for subsequent use in figure
% computes 
% - saccade frequency
% - fixation stability used for analysis of behavioral data


% continous data, [SampleTime,x,y,pupil,Vel_X,Vel_Y,Ve,Acc_x,Acc_y,Acc]
SubList = {'A01','A02','A04','A05','A06','A07','A08','A09','A10','A11','A12','A13','A14','A15','A16','A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A32','A33','A34','A35','A36','A37'};

WHICH_EXP = 4;

% ------------------------------------------
% parameters
ARG =[];
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset

% ---------------------------------------------------------------------------------
% load behavioral data 
Nsub = length(SubList);
defdirsCK15;  
sname = sprintf('%s/ProcessedData_Exp%d.mat',Prepropath,WHICH_EXP);
x = load(sname);
Data_Beh = x.Data_all;
ARG.VarList = x.ARG.VarList;

Ind_Freq =  find(strcmp(ARG.VarList,'Freq'));
Ind_Ear =  find(strcmp(ARG.VarList,'Ear'));
Ind_SOA =  find(strcmp(ARG.VarList,'SOA'));
Ind_SOABin =  find(strcmp(ARG.VarList,'SOABin'));
Ind_RT =  find(strcmp(ARG.VarList,ARG.WHICHRT));
Ind_Resp = find(strcmp(ARG.VarList,'Resp'));
Ind_Thr = find(strcmp(ARG.VarList,'SNR'));
Ind_TrialId = find(strcmp(ARG.VarList,'TotalTrialNr'));

% for dual task
Ind_task =  find(strcmp(ARG.VarList,'Timechange'));
Ind_stimchange = find(strcmp(ARG.VarList,'Fixchange'));
Ind_taskresp = find(strcmp(ARG.VarList,'ChangeResp'));


% ---------------------------------------------------------------------------------
% for visualization of longer tiem course, we use only trials wiht RT > 0.9
% ~ 50% of trials

Ntrials=zeros(Nsub,2);
AllData = zeros(20,6,190);
for S=1:Nsub

  % get behavioral data 
  Beh = Data_Beh{S}(:,[Ind_RT,Ind_Resp,Ind_SOA,Ind_Ear,Ind_TrialId,Ind_SOA]);

  % get eye data
  snameE =  sprintf('%s/Exp%d/eye/EyePrepro_%s.mat',Datapath,WHICH_EXP,SubList{S});
  dataE = load(snameE);
  fs = round(1./median(diff(dataE.tax_EyedataDispOn)));

  Ntrials(S,1) = size(dataE.EyedataDispOn,1);

  % ----------------------------------------------------
  % continous eye data
  EyeData = double(dataE.EyedataDispOn);
  % resample to 100 hz.
  EyeData = resample(EyeData,100,fs,'Dimension',3);
  taxnew = [dataE.tax_EyedataDispOn(1):0.01 :dataE.tax_EyedataDispOn(end)];
  J = find( (taxnew>-0.4).*(taxnew<=1.5));

  % ----------------------------------------------------
  % sort Beh and Eye trials in same order to facilitate selection
  EyeId = dataE.EyeId;
  [ia,ib,ic] = intersect(EyeId,Beh(:,5));
  Beh = Beh(ic,:);
  EyeId = EyeId(ib,:);
  % sort both in ascending order
  [~,o] = sort(Beh(:,5));
  Beh = Beh(o,:);
  [~,o] = sort(EyeId);
  EyeId = EyeId(o,:);
  EyeData = EyeData(o,:,:);
  
  Saccades = dataE.SaccTrials(o);

  % reduce to good trials
  isgood = sum(isnan(EyeData(:,1,:)),3)==0;
  Ntrials(S,2) =  sum(isgood);
  EyeData = EyeData(isgood,:,:);
  Beh = Beh(isgood,:);
  Saccades = Saccades(isgood);
  EyeId = EyeId(isgood,:);

  % ---------------------------------------------------------------------------------
  % for visualization
  % select trials with late RT as early RTs result in bad eye data. only
  % for the figure !
  islate = find(Beh(:,1)>0.9);

  % prepare dependent variables
  Vlist =[2,3,4]; %x,y,pupil
  c=1;
  for v=1:3
    tmp = squeeze(EyeData(islate,Vlist(v),:));
    tmp = tmp(:,J);
    % standardize pupil again
    if v==3
      tmp = tmp-mean(tmp(:),'omitnan');
      tmp = tmp./std(tmp(:),'omitnan');
    end
    % average
    AllData(S,c,:) =  mean(tmp);
    c=c+1;
    if v<3
      % std
      AllData(S,c,:) =  std(tmp);
      c=c+1;
    end
  end

  
  % compute saccade frequency on same time axis
  for t=1:length(Saccades)
    % get time of each saccade
    %  Eye.sacc  [start SecTime, end SecTime, duration SecTime, Xstart, Ystart, Xend, Yend,  avgvel, peakvel, start SampleTime, end SampleTime, pVel,PAcc]
    for l=1:size(Saccades{t},1)
      ton = Saccades{t}(l,1);
      [~,ind] = min(abs(taxnew(J)-ton));
      AllData(S,c,ind) =     AllData(S,c,ind)+1;
    end
  end
  % normalize by available number of trials
  AllData(S,c,:) = AllData(S,c,:) /length(Saccades);


  % ---------------------------------------------------------------------------------
  % measure of fixation stability and mean pupil during stim prior to
  % target


  % pupil
  pupil = squeeze(EyeData(:,4,:));
  % standardize pupil again
  pupil = pupil-mean(pupil(:),'omitnan');
  pupil = pupil./std(pupil(:),'omitnan');
  
  % generate saccade vector
  SaccadeV = zeros(length(Saccades),length(taxnew));
  for t=1:length(Saccades)
    % get time of each saccade
    %  Eye.sacc  [start SecTime, end SecTime, duration SecTime, Xstart, Ystart, Xend, Yend,  avgvel, peakvel, start SampleTime, end SampleTime, pVel,PAcc]
    for l=1:size(Saccades{t},1)
      ton = Saccades{t}(l,1);
      [~,ind] = min(abs(taxnew(J)-ton));
      SaccadeV(t,ind) = 1;
    end
  end


  [~,t0] = min(abs(taxnew-0));
  PerTrial=[];
  for t=1:length(Saccades)
    t_target = Beh(t,6)+0.3; % time of target rel to 0
    [~,t_target] = min(abs(taxnew-t_target));
    % pupil
    PerTrial(t,1) =  mean(pupil(t,[t0:t_target]));
    
    % fixation stability - mean of SDs of eye position 
    PerTrial(t,2) = mean(std(EyeData(t,[2,3],[t0:t_target]),[],3),2);

    % fixation stability - mean velocity
    PerTrial(t,3)  = mean(EyeData(t,7,[t0:t_target]),3);
     % fixation stability - mean velocity
    PerTrial(t,4)  = mean(SaccadeV(t,[t0:t_target]),2);
  end

  % for sanity check, store correlations with SOA
  CorrsPerTrials(S,1) = corr(PerTrial(:,1),Beh(:,6));
  CorrsPerTrials(S,2) = corr(PerTrial(:,2),Beh(:,6));
  CorrsPerTrials(S,3) = corr(PerTrial(:,3),Beh(:,6));
  CorrsPerTrials(S,4) = corr(PerTrial(:,4),Beh(:,6));

  % implement a median split for pupil and fixation stability 
  % and save trial indices. We need the original trial numbers to match with behavior data !!!
  % Beh(:,5) or EyeId
  for f=1:2
    m = median(PerTrial(:,f));
    SplitTrials{S}{f}{1} = EyeId(find(PerTrial(:,f)>m));
    SplitTrials{S}{f}{2} = EyeId(find(PerTrial(:,f)<m));
  end
end

% we include only participants with > 600 trials in total
isgood = find(Ntrials(:,2)>600); % index in original SubList

AllData = AllData(isgood,:,:);
SubListUse = SubList(isgood); % those sub Ids with good eye tracking data
SplitTrials = SplitTrials(isgood);
% save split trials for separate analysis
sname = sprintf('%s/Split_trials_Exp4.mat',Prepropath);
save(sname,'SplitTrials','SubListUse','isgood')

taxf = taxnew(J);

sname = sprintf('%s/Data_for_eyefigure.mat',Prepropath);
save(sname,'taxf','AllData')



%% Figure with eye data for those trials wiht target > 0.9 sec


tstr = {'Mean Xpos','Std Xpos','Mean Ypos','Std Ypos','Pupil','Saccade'}
figure('Position',[   94         233        1190         560])
vpos = [1,4,2,5,3,6];
ranges(1,:) = [-0.5 0.5];
ranges(2,:) = [0 5];
ranges(3,:) = [-0.5 0.5];
ranges(4,:) = [0 5];
ranges(5,:) = [-0.2 0.3];
ranges(6,:) = [0 0.08];

ylab = {'deg','deg','deg','deg','Z','per Trial'}
for v=1:6
  subplot(2,3,vpos(v));
  tmp = sq(AllData(:,v,:));

  % time 50 is 0ms
  if sum(v==[1,3,5])
    tmp = tmp-tmp(:,50);
  end
  plot(taxf,tmp,'k');
  hold on
  plot(taxf,mean(tmp),'color',[0.5,0.85,0.28],'LineWidth',2.5)
  xlim([-0.4 1]); set(gca,'XTick',[-0.2:0.2:1],'XTickLabel',[]);
  title(tstr{v})
  ylim([ranges(v,1) ranges(v,2)])
  if vpos(v)>3
    xlabel('Time [s]');
    set(gca,'XTick',[-0.2:0.2:1],'XTickLabel',[-0.2:0.2:1]);
  end
  ylabel(ylab{v})
end


ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',11);


snamef = sprintf('%s/Eye properties.jpg',Figurepath);
% print('-djpeg','-r600',snamef);

