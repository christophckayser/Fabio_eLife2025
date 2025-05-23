

% Prevalence of rythmicity in participant samples for experiments 1-3
% figure 6 and figure 7 (collapsed acrosx experiments) 

clear;

defdirsCK15

% significance levels for each method
% 
Methods = {'AIC Bins','VS Bins','Spectra','AIC Trials','VS Trials'};
PCRIT005 = [NaN    2.4500    2.1000  -15.0000    2.5000];
PCRIT001 = [NaN    3.2500    2.6000  -10.0000    3.2500];
% note that AIC bins and AIC trials are not considered in the paper. 

% the specific values were obtained from the simulations to result in a
% calibrated false positive rate.

% -----------------------------------------------------------------------
% figure properties
figure('Position',[   174         160        1800         900])
nRows = 3 ; nCols = 8 ;
% - Create grid of axes.
[blx, bly] = meshgrid( 0.08:0.9/nCols:0.92, 0.12:0.8/nRows:0.8 ) ;
bly = bly([nRows:-1:1],:);

% adjust to make figure nice
blx(:,[2]) = blx(:,[2])-0.014;
blx(:,[3]) = blx(:,[3])-0.022;
blx(:,[4]) = blx(:,[4])+0.01;
blx(:,[6]) = blx(:,[6])-0.01;

blx(:,[7]) = blx(:,[7])+0.015;
blx(:,[8]) = blx(:,[8])+0.005;


hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, 0.66/nCols, 0.66/nRows] ), blx, bly, 'UniformOutput', false ) ;

% -----------------------------------------------------------------
% part 1 - spectral results
% -----------------------------------------------------------------

% color per condition
Colorvector(1,:) = [40,40,40]/250; % both
Colorvector(2,:) = [250,150,40]/250; % left
Colorvector(3,:) = [40,150,250]/250; % right
valnames = {'d-prime','bias','rt'};
earname = {'Both','Left','Right'};

for WHICH_EXP = 1:3
  defdirsCK15;
  sname = sprintf('%s/Spectra_Exp%d.mat',Prepropath,WHICH_EXP);
  load(sname,'ARG','SpectraShuf','fax');

  sname = sprintf('%s/SpectraMC_Exp%d.mat',Prepropath,WHICH_EXP);
  load(sname,'SpectraMC');
  nsim = size(SpectraMC{1,1},1);

  % p-values for spectra
  pcrit005 = PCRIT005(3);
  pcrit001 = PCRIT001(3);
  for par=1:3
    axes( hAxes{WHICH_EXP,par} ); hold on
    for E=1:3
      SpecShuf  = squeeze(mean(SpectraShuf{par,E},1));
      
      % using the threshold we count the fraction of MC samples that reach
      % threshold for each frequency

      % for each frequency, find the VS in the shuffled data corresponding
      % to the respective percentile given by pcrit005
      ptile = 10^(-pcrit005);
      ptile = 100-ptile*100;
      thr = prctile(SpecShuf,ptile);
      
      y = 100*mean(SpectraMC{par,E}>repmat(thr,[nsim,1]));
      PrevalenceSpec{WHICH_EXP}(par,E,:) = y;
      plot(fax,y,'-','color',Colorvector(E,:),'LineWidth',1.2);
    end
    ylim([0 100]);  set(gca,'YTick',[0:20:100]);   set(gca,'YTickLabel',[]);
    xlim([1 8]);   set(gca,'XTick',[1:1:8]);    set(gca,'XTickLabel',[]);
    l = line([fax(1) fax(end)],[50 50 ],'LineStyle','--','color',ones(1,3)*0.5,'LineWidth',1);
    if par==1
      ylabel('% sig. samples')
       set(gca,'YTickLabel',[0:20:100]);
    end
    if WHICH_EXP==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[1:1:8]);
    elseif WHICH_EXP==1
      title(sprintf('%s',valnames{par}))
      if par==1
        text(1.4,90,earname{1},'color',Colorvector(1,:),'FontWeight','bold','FontSize',12)
        text(3.6,90,earname{2},'color',Colorvector(2,:),'FontWeight','bold','FontSize',12)
        text(5.7,90,earname{3},'color',Colorvector(3,:),'FontWeight','bold','FontSize',12)
      end
    end
  end
end
 
clear SpectraShuf SpectraMC SpecShuf ARG
fax_spec = fax;

% -----------------------------------------------------------------
% part 2 - binned statistics
% -----------------------------------------------------------------
% color per condition
Colorvector(1,:) = [40,40,40]/250; % both
Colorvector(2,:) = [250,150,40]/250; % left
Colorvector(3,:) = [40,150,250]/250; % right
valnames = {'d-prime','bias','rt'};

for WHICH_EXP = 1:3
  sname = sprintf('%s/Binned_Exp%d.mat',Prepropath,WHICH_EXP);
  X = load(sname,'VectorShuf','ARG');
  sname = sprintf('%s/BinnedMC_Exp%d.mat',Prepropath,WHICH_EXP);
  load(sname,'VectorMC')
  nsim = size(VectorMC{1,1},1);

  ARG = X.ARG;
  % p-value for vector strength of rhythmic component
  pcrit005 = PCRIT005(2);
  pcrit001 = PCRIT001(2);

  for par=1:3
    axes( hAxes{WHICH_EXP,par+3} ); hold on
    for E=1:3
      % vector length and group mean
      VectorShuf = squeeze(mean(X.VectorShuf{par,E},1));
      VectorShuf =VectorShuf(:,[2:end]);
      % using the threshold we count the fraction of MC samples that reach
      % threshold for each frequency

      % for each frequency, find the VS in the shuffled data corresponding
      % to the respective percentile given by pcrit005
      ptile = 10^(-pcrit005);
      ptile = 100-ptile*100;
      thr = prctile(VectorShuf,ptile);

      y = 100*mean(VectorMC{par,E}(:,[2:end])>repmat(thr,[nsim,1]));
      plot([ARG.flist ],y,'color',Colorvector(E,:),'LineWidth',1.2);
      PrevalenceBin{WHICH_EXP}(par,E,:) = y;
    end
    ylim([0 100]);  set(gca,'YTick',[0:20:100]);   set(gca,'YTickLabel',[]);
    xlim([1 8]);   set(gca,'XTick',[1:8]);    set(gca,'XTickLabel',[]);
    l = line([ARG.flist(1) ARG.flist(end)],[50 50 ],'LineStyle','--','color',ones(1,3)*0.5,'LineWidth',1);
    if par==1
      set(gca,'YTickLabel',[0:20:100]);
    end
    if WHICH_EXP==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[1:8]);
    elseif WHICH_EXP==1
      title(sprintf('%s',valnames{par}))
    end

  end
  drawnow
end
fax_bin = ARG.flist;

clear X  VectorMC VectorShuf ARG


%% -----------------------------------------------------------------
% part 3 - single trial statistics
% -----------------------------------------------------------------

valnames = {'accuracy','rt'};
% color per condition
Colorvector(1,:) = [250,150,40]/250; % left
Colorvector(2,:) = [40,150,250]/250; % right
Colorvector(3,:) = [40,40,40]/250; % both
earname = {'Left','Right','Both'};
PvalBeta=[];

for WHICH_EXP = 1:3
  sname = sprintf('%s/Singletrials_Exp%d.mat',Prepropath,WHICH_EXP);
  X = load(sname,'VectorShuf','ARG');
  ARG = X.ARG;
  sname = sprintf('%s/SingletrialsMC_Exp%d.mat',Prepropath,WHICH_EXP);
  load(sname,'VectorMC')
  nsim = size(VectorMC{1,1},1);

  pcrit005 = PCRIT005(5);
  pcrit001 = PCRIT001(5);
  for par=1:2 % RT/ accuracy
    axes( hAxes{WHICH_EXP,par+6} ); hold on
    for E=1:3 % ear
      % vector length and group mean
      VectorShuf = squeeze(mean(X.VectorShuf{par,E},1));
      VectorShuf = VectorShuf(:,[2:end]);
      % using the threshold we count the fraction of MC samples that reach
      % threshold for each frequency

      % for each frequency, find the VS in the shuffled data corresponding
      % to the respective percentile given by pcrit005
      ptile = 10^(-pcrit005);
      ptile = 100-ptile*100;
      thr = prctile(VectorShuf,ptile);

      y = 100*mean(VectorMC{par,E}(:,[2:end])>repmat(thr,[nsim,1]));
      PrevalenceTrials{WHICH_EXP}(par,E,:) = y;

      plot([ARG.flist ],y,'color',Colorvector(E,:),'LineWidth',1.2);
    end
    l = line([ARG.flist(1) ARG.flist(end)],[50 50 ],'LineStyle','--','color',ones(1,3)*0.5,'LineWidth',1);
    ylim([0 100]);  set(gca,'YTick',[0:20:100]);   set(gca,'YTickLabel',[]);
    xlim([1 12]);   set(gca,'XTick',[2:2:12]);    set(gca,'XTickLabel',[]);
    if WHICH_EXP==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[2:2:12])
    elseif WHICH_EXP==1
      title(sprintf('%s',valnames{par}))
    end
    if par==1
      set(gca,'YTickLabel',[0:20:100]);
    end
    xtickangle(0)
  end

end

clear X  VectorMC VectorShuf


ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',12);

% -----------------------------------------
% add experimentl labels
for WHICH_EXP=1:3
  axes( hAxes{WHICH_EXP,1} );
  text(-2.2,40,sprintf('Exp %d',WHICH_EXP),'FontSize',16,'Rotation',90);
end

% -----------------------------------------
% add analysis labels
axes( hAxes{1,2} );
text(4,120,'Spectra','FontSize',16);

axes( hAxes{1,5} );
text(4,120,'Binned','FontSize',16);

axes( hAxes{1,7} );
text(11,120,'Trials','FontSize',16);


% -----------------------------------------
% add panels
axes( hAxes{1,1} );
text(-1,120,'A','FontSize',20,'FontWeight','bold');

axes( hAxes{1,4} );
text(-1,120,'B','FontSize',20,'FontWeight','bold');

axes( hAxes{1,7} );
text(-1,120,'C','FontSize',20,'FontWeight','bold');

% -----------------------------------------
% add highlight for clear effects
pans = [1,4,7];
WHICH_EXP = 1;
for p=1:length(pans)
  axes( hAxes{WHICH_EXP,pans(p)} );
  h = fill([1.4 2.5 2.5 1.4],[0 0 100 100],[120 120 240]/240,'LineStyle','none');
  set(h,'FaceAlpha',0.3);
end

WHICH_EXP = 3;
pans = [1,4];
for p=1:length(pans)
  axes( hAxes{WHICH_EXP,pans(p)} );
  h = fill([2.3 3.4 3.4 2.3],[0 0 100 100],[120 120 120]/240,'LineStyle','none');
  set(h,'FaceAlpha',0.3);
end

WHICH_EXP = 2;
pans = [2,5,7];
for p=1:length(pans)
  axes( hAxes{WHICH_EXP,pans(p)} );
  h = fill([6 7 7 6],[0 0 100 100],[240 120 100]/240,'LineStyle','none');
  set(h,'FaceAlpha',0.3);
end

snamef = sprintf('%s/Figure6_prevalence.jpg',Figurepath);
print('-djpeg','-r600',snamef);







%% ----------------------------------------------------------
% consensus maps - Figure 6


Colorvector(1,:) = [40,40,40]/250; % both
Colorvector(2,:) = [250,150,40]/250; % left
Colorvector(3,:) = [40,150,250]/250; % right
valnames = {'d-prime','bias','rt'};
earname = {'Both','Left','Right'};

x = PrevalenceSpec{1};
for e=2:3
  x = x+PrevalenceSpec{e};
end
x = x/3;

figure('Position',[    102         207        1546         342
]);
nCols = 4 ;
[blx, bly] = meshgrid( 0.08:0.9/nCols:0.92, 0.18) ;
% adjust to make figure nice
blx(:,[2]) = blx(:,[2])-0.014;
blx(:,[3]) = blx(:,[3])+0.022;
%blx(:,[4]) = blx(:,[4])+0.01;
%blx(:,[6]) = blx(:,[6])-0.01;

% blx(:,[7]) = blx(:,[7])+0.015;
% blx(:,[8]) = blx(:,[8])+0.005;

LW = [0.8,1.3,1.3];
LS{1} = '--'; LS{2} ='-'; LS{3} ='-';
hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, 0.66/nCols, 0.65] ), blx, bly, 'UniformOutput', false ) ;

for par=1:2
  axes(hAxes{par}); hold on
  for E=1:3
    plot(fax_spec,sq(x(par,E,:)),LS{E},'color',Colorvector(E,:),'LineWidth',LW(E));
  end
  ylim([0 40]);  set(gca,'YTick',[0:10:40]);   set(gca,'YTickLabel',[]);
  xlim([1 8]);   set(gca,'XTick',[1:1:8]);    set(gca,'XTickLabel',[]);
  title(sprintf('%s',valnames{par}));
  xlabel('Frequency [Hz]');
  set(gca,'XTickLabel',[1:8]);

  if par==1
    set(gca,'YTickLabel',[0:10:40]);;
    ylabel('# sig. samples');
    text(1.4,37,earname{1},'color',Colorvector(1,:),'FontWeight','bold','FontSize',10)
    text(3.6,37,earname{2},'color',Colorvector(2,:),'FontWeight','bold','FontSize',10)
    text(5.7,37,earname{3},'color',Colorvector(3,:),'FontWeight','bold','FontSize',10)
  end
end

x = PrevalenceBin{1};
for e=2:3
  x = x+PrevalenceBin{e};
end
x = x/3;

for par=1:2
    axes(hAxes{par+2}); hold on

  for E=1:3
    plot(fax_bin,sq(x(par,E,:)),LS{E},'color',Colorvector(E,:),'LineWidth',LW(E));
  end
  ylim([0 40]);  set(gca,'YTick',[0:10:40]);   set(gca,'YTickLabel',[]);
  xlim([1 8]);   set(gca,'XTick',[1:1:8]);    set(gca,'XTickLabel',[]);
  title(sprintf('%s',valnames{par}));
  xlabel('Frequency [Hz]');
  set(gca,'XTickLabel',[1:8]);
  if par==1
    set(gca,'YTickLabel',[0:10:40]);;
    ylabel('# sig. samples');
  end
end



% -----------------------------------------
% add analysis labels
axes(hAxes{1});
text(8,47,'Spectra','FontSize',16);
text(-1,47,'A','FontSize',20,'FontWeight','bold');

axes(hAxes{3});
text(8,47,'Binned','FontSize',16);
text(-1,47,'B','FontSize',20,'FontWeight','bold');

% axes(hAxes{7});
% text(11,47,'Trials','FontSize',16);
% text(-1,47,'C','FontSize',20,'FontWeight','bold');
% 

ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',11);

snamef = sprintf('%s/Figure7_PrevalanceComb.jpg',Figurepath);
print('-djpeg','-r600',snamef);


