

clear;

% analyze task-thresholds 

defdirsCK15;

%  Sound ~ Tone_scale
%  Tone_scale = 10.^((Tone_SNR-4));
%  Tone_SNR   = THR.thresholds;

% -> Sound ~  10^THR * 10-4

% relative thresholds in dB is computed as 20 * log10(V/V0)
% -> 20*log10 (  10^THR * 10-4 / 10^THR_base * 10-4)
%    20*log10 ( 10^(Thr-Thr_base)) 
%    20*(Thr-Thr_base)


% note that this is voltage level, not presented or perceived sound!


figure(1);clf
for WHICH_EXP = 1:4

  sname = sprintf('%s/ProcessedData_Exp%d.mat',Prepropath,WHICH_EXP);
  x = load(sname);
  ARG.VarList = x.ARG.VarList;

  THR_stat=[];
  THR_time=[];

  for S=1:length(x.Data_all)
    tmp = x.Data_all{S};

    for f=1:2
      for e=1:2
        j = find( (tmp(:,2)==f).*(tmp(:,3)==e));
        % all trials for this (f,e)
        data = 20*tmp(j,9);
        
        % compute min/max, start , end 
        % because the actual value is somewhat arbirtary we normalize to
        % start value 
        offset = data(1);
        THR_stat(S,f,e,:) = [data(end)-offset,max(data)-min(data)];
        % and time course in 10 bins
        n = size(data,1);
        nbin = floor(n/10);
        data = data([1:nbin*10]);
        data = reshape(data,[nbin,10]);
        tmp2 = mean(data,1);
        THR_time(S,f,e,:) = tmp2-offset;

      end
    end
  end

  subplot(2,2,WHICH_EXP)
  %plot(sq(mean(mean(THR_time,2),3))','k')

  plot(sq(THR_time(:,1,2,:))','k')

  fprintf('Exp %d \n',WHICH_EXP)
  % statisstics on change 
  % change from start to end. averaged over frequencies and ears
  delta = mean(mean(THR_stat(:,:,:,1),2),3);
  % paired t-test
  [h,p,ci,stats] = ttest(delta,0);
  fprintf('start to end %1.2f+%1.2f [dB] \t p=%1.3f  t=%1.3f \n',mean(delta),sem(delta),p,stats.tstat)


end
