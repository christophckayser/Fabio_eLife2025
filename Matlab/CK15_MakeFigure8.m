

clear

% Results For exp 4, display a comparison of effect measures for
% - dual task
% - data split by pupil size 
% - data split by fixation stabily 
% we show the delta vector strength and the prevalence 

% SPLIT 1 : WITH DUAL TASK - WITHOUT
% SPLIT 2 : LARGE PUPIL - SMALL 
% SPLIT 3 : MORE EYE SD - LESS

% Figure 8

WHICH_EXP = 4;
defdirsCK15;

% color per condition
Colorvector(1,:) = [40,40,40]/250; % both
Colorvector(2,:) = [250,150,40]/250; % left
Colorvector(3,:) = [40,150,250]/250; % right
valnames = {'d-prime','bias','rt'};
earname = {'Both','Left','Right'};

USE_E = 1; % which ear to show. Here we show both.
% 1 is both, 2 is left, 3 is right

% significance lines. We run paired t-tests
% and correct for MCP across frequencies using FDR
pcrit = 0.05;


figure('Position',[    74         188        1093         810])
nRows = 3 ; nCols = 4 ;
% - Create grid of axes.
[blx, bly] = meshgrid( 0.18:0.8/nCols:0.9, 0.1:0.85/nRows:0.8 ) ;
bly = bly([nRows:-1:1],:);
blx(:,[2]) = blx(:,[2])-0.02;
blx(:,[3]) = blx(:,[3])+0.02;

hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, 0.66/nCols, 0.61/nRows] ), blx, bly, 'UniformOutput', false ) ;

% -----------------------------------------------------------------------
% comparison of delta vector strength
% -----------------------------------------------------------------------


% -----------------------------------------------------------------------
% DUAL TASK
sname = sprintf('%s/Timestats_Conditions_Exp%d.mat',Prepropath,WHICH_EXP);
X = load(sname);
ARG = X.ARG;
% compute condition differences
SPLIT=1;
for par=1:3 % parameter
  for E=1:3 % each ear, or combined
    DeltaVS{SPLIT,par,E} = X.Vector{1}{par,E}(:,[2:end])-X.Vector{2}{par,E}(:,[2:end]);
  end
end

% -----------------------------------------------------------------------
% EYE PROPERTIES 
sname = sprintf('%s/Timestats_eyesplit_Exp%d.mat',Prepropath,WHICH_EXP);
X = load(sname);

% compute condition differences
for SPLIT=1:2
  for par=1:3 % parameter
    for E= 1:3% each ear, or combined
      DeltaVS{SPLIT+1,par,E} = X.Vector{1,SPLIT}{par,E}(:,[2:end])-X.Vector{2,SPLIT}{par,E}(:,[2:end]);
    end
  end
end

Colorvector2 = [0.3,0.3,0.3];
Colorvector2 = repmat(Colorvector2,[10,1]);

for SPLIT=1:3
  for par=1:2
    % -----------------------------------------------------------------
    % vector strength
    axes( hAxes{SPLIT,par} ); hold on;

    for E=USE_E
      tmp =  DeltaVS{SPLIT,par,E};
      ckerrorshade(ARG.flist,mean(tmp),sem(tmp),'k');
      %  plot(ARG.flist,mean(tmp),'Color',Colorvector(E,:),'LineWidth',1.5);
      l = line([ARG.flist(1) ARG.flist(end)],[0 0],'color',ones(1,3)*0.5,'LineWidth',1);
      set(l,'HandleVisibility','off');
      df = size(tmp,1);
      pval = 2*tcdf(-abs(sqrt(df)*mean(tmp)./std(tmp)),df-1);
      [h, crit_p, adj_ci_cvrg, pval]=ck_stat_fdr(pval,0.05);
      issig = find(pval<pcrit);
      if ~isempty(issig)
        fprintf('VS %s %s  ',valnames{par},earname{E});
        for l=1:length(issig)
          fprintf('%2.1f   ',ARG.flist(issig(l)));
          plot(ARG.flist(issig(l)),mean(tmp(:,issig(l))),'.','MarkerSize',25,'Color','k');
        end
        fprintf('\n');
      end
    end
    xlim([0.5 8.5]); set(gca,'XTick',[2:2:12]);    set(gca,'XTickLabel',[])
    if SPLIT==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[2:2:12])
    end
    if par==1
      ylim([-0.15 0.15]);      set(gca,'YTick',[-0.1:0.1:0.1]);
      ylabel('\Delta Vector strength ');
    elseif par==2
      ylim([-0.15 0.15]);      set(gca,'YTick',[-0.1:0.1:0.1],'YTickLabel',[]);
    elseif par==3
      ylim([-0.015 0.015]);      set(gca,'YTick',[-0.01:0.01:0.01]);

    end
    if SPLIT==1
      title(sprintf('%s',valnames{par}))
    end
  end
end



% -----------------------------------------------------------------------
% prevalence
% -----------------------------------------------------------------------




% -----------------------------------------------------------------------
% DUAL TASK
sname = sprintf('%s/Timestats_ConditionsMC_Exp%d.mat',Prepropath,WHICH_EXP);
X = load(sname);
ARG = X.ARG;
% we have the uncorrected p-valus for each MC run
% apply correction for multiple tests
SPLIT=1;
for par=1:3 % parameter
  for E=1:3 % each ear, or combined
    nrun = size(X.PvalMC{par,E},1);
    Prev = zeros(nrun,size(X.PvalMC{par,E},2));
    for R=1:nrun
      [h, crit_p, adj_ci_cvrg, pval]=ck_stat_fdr(X.PvalMC{par,E}(R,:) ,0.05);
      Prev(R,:) = pval;
    end
    DeltaVS{SPLIT,par,E} = mean(Prev<pcrit)*100;
  end
end


% -----------------------------------------------------------------------
% EYE PROPERTIES 
sname = sprintf('%s/Timestats_eyesplitMC_Exp%d.mat',Prepropath,WHICH_EXP);
X = load(sname);

% compute condition differences
for SPLIT=1:2
  for par=1:3 % parameter
    for E= 1:3% each ear, or combined
      nrun = size(X.PvalMC{SPLIT}{par,E},1);
      Prev = zeros(nrun,size(X.PvalMC{SPLIT}{par,E},2));
      for R=1:nrun
        [h, crit_p, adj_ci_cvrg, pval]=ck_stat_fdr(X.PvalMC{SPLIT}{par,E}(R,:) ,0.05);
        Prev(R,:) = pval;
      end
      DeltaVS{SPLIT+1,par,E} = mean(Prev<pcrit)*100;
    end
  end
end

Colorvector2 = [0.3,0.3,0.3];
Colorvector2 = repmat(Colorvector2,[10,1]);
% ------------------------------------------------------------
% display
LW = [1.5,1,1];
for SPLIT=1:3
  for par=1:2
    axes( hAxes{SPLIT,par+2} ); hold on;
    for E=1:3
      plot(ARG.flist, DeltaVS{SPLIT,par,E},'color',Colorvector(E,:),'LineWidth',LW(E));
      l = line([ARG.flist(1) ARG.flist(end)],[50 50 ],'LineStyle','--','color',ones(1,3)*0.5,'LineWidth',1);
    end
    xlim([0.5 8.5]); set(gca,'XTick',[2:2:12]);    set(gca,'XTickLabel',[])
    ylim([0 100]); set(gca,'YTick',[0:20:100]);    set(gca,'YTickLabel',[])
    if SPLIT==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[2:2:12])
    end
    if par==1
      ylabel('% sig. samples');  set(gca,'YTickLabel',[0:20:100])
    end
    if SPLIT==1
      title(sprintf('%s',valnames{par}))
    end
  end
end







ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',12);


% add condition labels
tex = {'Dual task','Pupil size','Fixation stability'};
SPLIT=1;
axes( hAxes{SPLIT,2} );
text(7,0.25,sprintf('%s',tex{SPLIT}),'FontSize',16,'FontAngle','italic');

SPLIT=2;
axes( hAxes{SPLIT,2} );
text(7,0.2,sprintf('%s',tex{SPLIT}),'FontSize',16,'FontAngle','italic');

SPLIT=3;
axes( hAxes{SPLIT,2} );
text(6,0.2,sprintf('%s',tex{SPLIT}),'FontSize',16,'FontAngle','italic');


% -----------------------------------------
% add analysis labels
axes( hAxes{1,1} );
text(-2,0.18,'A','FontSize',20,'FontWeight','bold');

axes( hAxes{2,1} );
text(-2,0.18,'B','FontSize',20,'FontWeight','bold');

axes( hAxes{3,1} );
text(-2,0.18,'C','FontSize',20,'FontWeight','bold');


snamef = sprintf('%s/Figure8_Exp4.jpg',Figurepath);
print('-djpeg','-r600',snamef);

return


