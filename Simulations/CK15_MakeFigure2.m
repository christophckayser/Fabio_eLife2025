

clear;

% visualize results from simulations

% here we determine significance levels for each method.
% these have been determined below to adjust the FPR across simulations 
% with correction across frequencies to about 0.05 for each method


% resulting level of 0.05.
% these are just pasted here for visualization in the figure.
% they are actually computed below
PCRIT = [  NaN    2.4500    2.1000  -17.0000    2.5000];

%  Based on the code the methods are
Methods = {'AIC Binned','Binned','Spectra','Trials AIC','Trials'};


% with following colors
Colorvector(1,:) = [250,150,0]/250; %  AIC bins. Not shown. 
Colorvector(2,:) = [250,150,0]/250; % VS bins orange
Colorvector(3,:) = [80,80,80]/250; % Spectrum gray
Colorvector(4,:) = [50,150,200]/250; % AIC trials blue
Colorvector(5,:) = [0,200,80]/250; %  VS trials green
Titles = {'True positive rate','False positive rate','Asymmetrie','Variability'}

Nlev{1} = [2:1:8];
Nlev{2} = [2:2:6];
% -------------------------------------------------------------------
% binned data and spectra
% -------------------------------------------------------------------
clear UseSet*
fname = 'New2Simulations_bins.mat';
load(fname,'Set1','Set2','ARG');

%  GLM results:
%    out{1} AIC
%    out{2} PvalBeta;
%  Spectra results:
%    out{3} = PvalSpecAR;
%    out{4} = fax_spec;


% loop conditions
for UseSet=1:2
  Data=[];
  eval(sprintf('Data = Set%d',UseSet));
  nrun = size(Data{1},2);

  for n=1:length(Data)
    for rep=1:nrun
      aic =  Data{n}{rep}{1};
      % ---------------------------------------------------------------
      %  AIC difference for best frequency vs. null model
      aic = -(aic(2:end)-aic(1));

      Delta_AICBins{UseSet}(n,rep,:) = aic; % here positive AIC are indicatige of rhythmic effect

      % ---------------------------------------------------------------
      % p value vector strength
      Delta_VSBins{UseSet}(n,rep,:) =  -log10(Data{n}{rep}{2}+eps);
     
      % ---------------------------------------------------------------
      % p values for frequencies in spectrum
      Delta_AR{UseSet}(n,rep,:) = -log10(Data{n}{rep}{3}+eps);

    end
  end
end
fax_spec = Data{1}{1}{4};
fax_bins = ARG.flist;

% -------------------------------------------------------------------
% trial-based analysis
% -------------------------------------------------------------------

clear UseSet* ARG
fname = 'New2Simulations_trials1.mat';
load(fname,'Set1','ARG');
fname = 'New2Simulations_trials2.mat';
load(fname,'Set2');

for UseSet=1:2
  Data=[];
  eval(sprintf('Data = Set%d',UseSet));
  nrun = size(Data{1},2);

  for n=1:length(Data)
    for rep=1:nrun
      aic =  Data{n}{rep}{1};
      % ---------------------------------------------------------------
      %  AIC difference for best frequency vs. null model
      aic = -(aic(2:end)-aic(1));
      Delta_AICTrials{UseSet}(n,rep,:) = aic; % here positive AIC are indicatige of rhythmic effect

      % ---------------------------------------------------------------
      % p value vector strength
      Delta_VSTrials{UseSet}(n,rep,:) =  -log10(Data{n}{rep}{2}+eps);
    end
  end
end
fax_trials = ARG.flist;



%% -------------------------------------------------------
% rate of detecting effects at any frequency
% -------------------------------------------------------

figure(3);clf;

for UseSet=1:2
  subplot(2,2,UseSet); hold on

  % find signifciant effect across frequency per run

  plot(Nlev{UseSet},mean(sum( Delta_VSBins{UseSet}>PCRIT(2),3)>1,2),'color',Colorvector(2,:),'LineWidth',1.8);
  plot(Nlev{UseSet},mean(sum( Delta_AR{UseSet}>PCRIT(3),3)>1,2),'color',Colorvector(3,:),'LineWidth',1.8);
  plot(Nlev{UseSet},mean(sum( Delta_VSTrials{UseSet}>PCRIT(5),3)>1,2),'color',Colorvector(5,:),'LineWidth',1.8);

  ylabel('Significant findings')
  title(Titles{UseSet});

  if UseSet==1
    xlim([1.5 8.5]); ylim([-0.05 1.05]);
    set(gca,'XTick',[2:8]); set(gca,'YTick',[0:0.2:1]);
  elseif UseSet==2
    xlim([1.5 6.5]); ylim([0 0.2]);  set(gca,'YTick',[0:0.05:0.2]);
    set(gca,'XTick',[2,4,6]);
  end
  xlabel('SNR')
  grid on

  if UseSet==1
    posme = [0.1,0,0.4,0.2];
    for m=[2,3,5]
      text(2,0.8-posme(m-1),Methods{m},'color',Colorvector(m,:),'FontWeight','Bold','FontSize',12);
    end
  end
end

% -------------------------------------------------------
% display 'significance' vs. frequency
UseSet=1;

subplot(2,3,5);
plot(fax_bins,sq(mean(Delta_VSBins{UseSet},2))','color',Colorvector(2,:));
xlabel('frequency');
ylabel('-log10(p)')
xlim([1 8])
title(Methods{2})

subplot(2,3,4);
plot(fax_spec,sq(mean(Delta_AR{UseSet},2))','color',Colorvector(3,:));
xlabel('frequency');
ylabel('-log10(p)')
xlim([1 8])
title(Methods{3})

% subplot(2,4,8);
% plot(fax_trials,sq(mean(Delta_AICTrials{UseSet},2))','color',Colorvector(4,:));
% ylabel('\Delta-AIC')
% xlabel('frequency');
% xlim([1 12])
% title(Methods{4})


subplot(2,3,6);
plot(fax_trials,sq(mean(Delta_VSTrials{UseSet},2))','color',Colorvector(5,:));
ylabel('-log10(p)')
xlabel('frequency');
xlim([1 12])
title(Methods{5})





% -----------------------------------------
% add  labels
subplot(2,2,1);
text(0.45,1.15,'A','FontSize',20,'FontWeight','bold');

subplot(2,2,2); hold on
text(0.3,0.22,'B','FontSize',20,'FontWeight','bold');


subplot(2,2,2);
pos = get(gca,'Position')
pos(3) = pos(3)-0.1;
set(gca,'Position',pos)

subplot(2,3,4); hold on
text(-1.2,16,'C','FontSize',20,'FontWeight','bold');




set(gcf,'Position',[  205 214 1050 750])


ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',11);

Figurepath='F:/CKDATA/Projects/projects/Hearing/CK15/figures';

snamef = sprintf('%s/Figure2_simulations.jpg',Figurepath);
print('-djpeg','-r600',snamef);



% -------------------------------------------------------------------------------
% compute FPR for different thresholds to obtain the final thresholds to be
% used
% to choose same actual FPR across methods
% -------------------------------------------------------------------------------

% Methods = {'AIC Bins','VS Bins','Spectra','AIC Trials','VS Trials'};

% testing levels which we probe to compute the FPR
P2list = [1:0.05:3.5]; % first level p values
P1list = [-30:0.2:2]; % first level delta AIC  values

for p=1:length(P2list)
  % average across noise levels
  fpr(p,2) =  mean(mean(sum( Delta_VSBins{2}>P2list(p),3)>1,2));
  fpr(p,3) =  mean(mean(sum( Delta_AR{2}>P2list(p),3)>1,2));
  fpr(p,5) =  mean(mean(sum( Delta_VSTrials{2}>P2list(p),3)>1,2));
end
% for AIC
for p=1:length(P1list)
  fpraic(p) =  mean(mean(sum( Delta_AICTrials{2}>P1list(p),3)>1,2));
end
% to find the threshold yielding a FPR of 0.05
PTARGET = 0.05;
[~,o] = min(abs(fpr-PTARGET),[],1);
[~,o2] = min(abs(fpraic'-PTARGET),[],1);
PCRIT005 = [NaN,P2list(o(2)),P2list(o(3)), P1list(o2),P2list(o(5))]
% to find the threshold yielding a FPR of 0.01
PTARGET = 0.01;
[~,o] = min(abs(fpr-PTARGET),[],1);
[~,o2] = min(abs(fpraic'-PTARGET),[],1);

PCRIT001 = [NaN,P2list(o(2)),P2list(o(3)), P1list(o2),P2list(o(5))]


% compute sensitivity and specificit for this choice then


%Sensitivity (true positive rate) is the probability of a positive test result, conditioned on the individual truly being positive.
%Specificity (true negative rate) is the probability of a negative test result, conditioned on the individual truly being negative.

Sensitivity{2} = mean(sum( Delta_VSBins{1}>PCRIT(2),3)>1,2);
Sensitivity{3} = mean(sum( Delta_AR{1}>PCRIT(3),3)>1,2);
Sensitivity{4} = mean(sum( Delta_AICTrials{1}>PCRIT(4),3)>1,2);
Sensitivity{5} = mean(sum( Delta_VSTrials{1}>PCRIT(5),3)>1,2);
% reduce to matching SNRs
for me =2:5
  Sensitivity{me} = Sensitivity{me}([1,2,4]);
end

Specificity{2} = 1-mean(sum( Delta_VSBins{2}>PCRIT(2),3)>1,2);
Specificity{5} = 1-mean(sum( Delta_AR{2}>PCRIT(3),3)>1,2);
Specificity{3} = 1-mean(sum( Delta_AICTrials{2}>PCRIT(4),3)>1,2);
Specificity{4} = 1-mean(sum( Delta_VSTrials{2}>PCRIT(5),3)>1,2);

% average noise levels
for me=2:5
  Outcome(me,:) = [mean(Sensitivity{me}),mean(Specificity{me})];
end
Outcome



