
% show statistical results (p-values) for experiments 1-3.
% Figure 4 

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
% part 1 - spectral analysis
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
  load(sname,'ARG','SpectraShuf','Spectra','fax');

  pcrit005 = PCRIT005(3);
  pcrit001 = PCRIT001(3);
  for par=1:3
    axes( hAxes{WHICH_EXP,par} ); hold on
    for E=1:3
      % compute group- averages for actual and surrogate data
      Spectrue = squeeze(mean(Spectra{par,E}(:,:),1));
      SpecShuf  = squeeze(mean(SpectraShuf{par,E},1));
      PvalSpecAR=[];
      % compute percentile p-values
      for f=1:length(fax)
        PvalSpecAR(f) = mean(SpecShuf(:,f)>Spectrue(f));
      end
      y = -log10(PvalSpecAR+eps);
      plot(fax,y,'-','color',Colorvector(E,:),'LineWidth',1.2);

      issig = find(y>pcrit005);
     if length(issig)>1
        fprintf('Exp %d %s %s  ',WHICH_EXP,valnames{par},earname{E});
        for l=1:length(issig)
          fprintf('%2.1f   ',fax(issig(l)));
          plot(fax(issig(l)),y(issig(l)),'.','color',Colorvector(E,:),'MarkerSize',18);
        end
        fprintf('\n');
      end
    end
    l = line([fax(1) fax(end)],[pcrit005 pcrit005],'color',ones(1,3)*0.5,'LineWidth',1.2);
    l = line([fax(1) fax(end)],[pcrit001 pcrit001],'color',ones(1,3)*0.5,'LineWidth',0.8,'LineStyle','--');
    set(l,'HandleVisibility','off');
    ylim([0 4]);  set(gca,'YTick',[0:1:4]);   set(gca,'YTickLabel',[]);
    xlim([1 8]);   set(gca,'XTick',[1:1:8]);    set(gca,'XTickLabel',[]);
    if par==1
      ylabel('-log10(p)')
      set(gca,'YTickLabel',[0:4]);
    end
    if WHICH_EXP==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[1:1:8]);

    elseif WHICH_EXP==1
      title(sprintf('%s',valnames{par}))
      if par==1
        text(1.4,3.6,earname{1},'color',Colorvector(1,:),'FontWeight','bold','FontSize',12)
        text(3.6,3.6,earname{2},'color',Colorvector(2,:),'FontWeight','bold','FontSize',12)
        text(5.7,3.6,earname{3},'color',Colorvector(3,:),'FontWeight','bold','FontSize',12)
      end
    end
  end
end % EXP

clear SpectraShuf Spectra SpecShuf ARG

% -----------------------------------------------------------------
% part 2 - binned statistics
% color per condition
Colorvector(1,:) = [40,40,40]/250; % both
Colorvector(2,:) = [250,150,40]/250; % left
Colorvector(3,:) = [40,150,250]/250; % right
valnames = {'d-prime','bias','rt'};

for WHICH_EXP = 1:3
  sname = sprintf('%s/Binned_Exp%d.mat',Prepropath,WHICH_EXP);
  X = load(sname,'Fit_Qual','Vector','VectorShuf','ARG');
  ARG = X.ARG;
  % p-value for vector strength of rhythmic component
  pcrit005 = PCRIT005(2);
  pcrit001 = PCRIT001(2);

  for par=1:3
    axes( hAxes{WHICH_EXP,par+3} ); hold on
    for E=1:3
      % compute group- averages for actual and surrogate data
      Vector = mean(X.Vector{par,E},1);
      VectorShuf = squeeze(mean(X.VectorShuf{par,E},1));
      % consider only frequency dependent models
      Vector = Vector([2:end]);
      VectorShuf = VectorShuf(:,[2:end]);
      PvalBeta=[];
      % compute percentile p-values
      for f=1:length(ARG.flist)
        PvalBeta(f) = mean(VectorShuf(:,f)>Vector(f));
      end
      y = -log10(PvalBeta+eps)';
      plot([ARG.flist ],y,'color',Colorvector(E,:),'LineWidth',1.2);
      issig = find(y>pcrit005);
      if length(issig)>1
        fprintf('Exp %d %s %s  ',WHICH_EXP,valnames{par},earname{E});
        for l=1:length(issig)
          fprintf('%2.1f   ',ARG.flist(issig(l)));
          plot(ARG.flist(issig(l)),y(issig(l)),'.','color',Colorvector(E,:),'MarkerSize',18);
        end
        fprintf('\n');
      end
    end
    l = line([fax(1) fax(end)],[pcrit005 pcrit005],'color',ones(1,3)*0.5,'LineWidth',1.2);
    l = line([fax(1) fax(end)],[pcrit001 pcrit001],'color',ones(1,3)*0.5,'LineWidth',0.8,'LineStyle','--');
    set(l,'HandleVisibility','off');
    ylim([0 4]);  set(gca,'YTick',[0:1:4]);   set(gca,'YTickLabel',[]);
    xlim([1 8]);   set(gca,'XTick',[1:8]);    set(gca,'XTickLabel',[]);
   
    if par==1
      set(gca,'YTickLabel',[0:4]);
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

clear Fit_Qual Vector VectorShuf X ARG

% -----------------------------------------------------------------
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
  X = load(sname,'Vector','VectorShuf','ARG');
  ARG = X.ARG;
  % ------------------------------------------------------
  % p-value for vector strength of rhythmic component
  % vector length and group mean
  for par=1:2 % RT/ accuracy
    for E=1:3 % ear
      % compute group- averages for actual and surrogate data

      Vector = squeeze(mean(X.Vector{par,E},1))';
      VectorShuf = squeeze(mean(X.VectorShuf{par,E},1));
      % consider only frequency dependent models

      Vector = Vector([2:end]);
      VectorShuf = VectorShuf(:,[2:end]);
      PvalBeta{par,E} =[];
      % compute percentile p-values

      for f=1:length(ARG.flist)
        PvalBeta{par,E}(f) = mean(VectorShuf(:,f)>Vector(f));
      end

    end
  end
  % p-value for vector strength of rhythmic component
  pcrit005 = PCRIT005(5);
  pcrit001 = PCRIT001(5);
  for par=1:2 % RT/ accuracy
     axes( hAxes{WHICH_EXP,par+6} ); hold on
    for E=1:3 % ear
      y = -log10(PvalBeta{par,E}+eps)';
      plot([ARG.flist ],y,'color',Colorvector(E,:),'LineWidth',1.2);
      issig = find(y>pcrit005);
       if length(issig)>1
        fprintf('Exp %d %s %s  ',WHICH_EXP,valnames{par},earname{E});
        for l=1:length(issig)
          fprintf('%2.1f   ',ARG.flist(issig(l)));
          plot(ARG.flist(issig(l)),y(issig(l)),'.','color',Colorvector(E,:),'MarkerSize',18);
        end
        fprintf('\n');
      end
    end
    l = line([ARG.flist(1) ARG.flist(end)],[pcrit005 pcrit005],'color',ones(1,3)*0.5,'LineWidth',1.2);
    l = line([ARG.flist(1) ARG.flist(end)],[pcrit001 pcrit001],'color',ones(1,3)*0.5,'LineWidth',0.8,'LineStyle','--');
    set(l,'HandleVisibility','off');
    ylim([0 4]);  set(gca,'YTick',[0:1:4]);   set(gca,'YTickLabel',[]);
    xlim([1 12]);   set(gca,'XTick',[2:2:12]);    set(gca,'XTickLabel',[]);
    if WHICH_EXP==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[2:2:12])
    elseif WHICH_EXP==1
      title(sprintf('%s',valnames{par}))
    end
    if par==1
      set(gca,'YTickLabel',[0:4]);
    end
    xtickangle(0)
  end

end

clear Fit_Qual Vector VectorShuf X


ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',12);

% -----------------------------------------
% add experimentl labels
for WHICH_EXP=1:3
  axes( hAxes{WHICH_EXP,1} );
  text(-2.2,1.6,sprintf('Exp %d',WHICH_EXP),'FontSize',16,'Rotation',90);
end

% -----------------------------------------
% add analysis labels
axes( hAxes{1,2} );
text(4,5.2,'Spectra','FontSize',16);

axes( hAxes{1,5} );
text(4,5.2,'Binned','FontSize',16);

axes( hAxes{1,7} );
text(11,5.2,'Trials','FontSize',16);


% -----------------------------------------
% add panels
axes( hAxes{1,1} );
text(-1,5,'A','FontSize',20,'FontWeight','bold');

axes( hAxes{1,4} );
text(-1,5,'B','FontSize',20,'FontWeight','bold');

axes( hAxes{1,7} );
text(-1,5,'C','FontSize',20,'FontWeight','bold');

% -----------------------------------------
% add highlight for clear effects
pans = [1,4,7];
WHICH_EXP = 1;
for p=1:length(pans)
  axes( hAxes{WHICH_EXP,pans(p)} );
  h = fill([1.4 2.5 2.5 1.4],[0 0 4 4],[120 120 240]/240,'LineStyle','none');
  set(h,'FaceAlpha',0.3);
end

WHICH_EXP = 3;
pans = [1,4];
for p=1:length(pans)
  axes( hAxes{WHICH_EXP,pans(p)} );
  h = fill([2.3 3.4 3.4 2.3],[0 0 4 4],[120 120 120]/240,'LineStyle','none');
  set(h,'FaceAlpha',0.3);
end

WHICH_EXP = 2;
pans = [2,5,7];
for p=1:length(pans)
  axes( hAxes{WHICH_EXP,pans(p)} );
  h = fill([6 7 7 6],[0 0 4 4],[240 120 100]/240,'LineStyle','none');
  set(h,'FaceAlpha',0.3);
end

snamef = sprintf('%s/Figure4_pvals.jpg',Figurepath);
print('-djpeg','-r600',snamef);


