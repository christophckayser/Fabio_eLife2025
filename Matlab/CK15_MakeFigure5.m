
% show the actual effects for the spectral analysis
% including the shuffled estimates.
% Figure 5


clear; close all

defdirsCK15

WHICH_EXP = 1; % only exp 1
% significance levels for each method
METHOD = 2;
Methods = {'AIC Bins','VS Bins','Spectra','AIC Trials','VS Trials'};
PCRIT005 = [NaN    2.4500    2.1000  -15.0000    2.5000];
PCRIT001 = [NaN    3.2500    2.6000  -10.0000    3.2500];
% note that AIC bins and AIC trials are not considered in the paper. 
% the specific values were obtained from the simulations to result in a
% calibrated false positive rate.


% -----------------------------------------------------------------------
% figure properties
figure('Position',[  150   196   969   803])
nRows = 3 ; nCols = 3 ;
% - Create grid of axes.
[blx, bly] = meshgrid( 0.14:0.9/nCols:0.9, 0.1:0.8/nRows:0.8 ) ;
bly = bly([nRows:-1:1],:);
%blx(:,[2]) = blx(:,[2])-0.07;
% blx(:,[3]) = blx(:,[3])-0.14;

ColSize = [0.19 0.19 0.19 ];
ColSize = repmat(ColSize,[3,1]);
hAxes = arrayfun( @(x,y,z) axes( 'Position', [x, y, z, 0.66/nRows] ), blx, bly, ColSize, 'UniformOutput', false ) ;

% color per condition
Colorvector(1,:) = [40,40,40]/250; % both
Colorvector(2,:) = [250,150,40]/250; % left
Colorvector(3,:) = [40,150,250]/250; % right
valnames = {'d-prime','bias','rt'};
earname = {'Both','Left','Right'};

Ymax =[0.4,0.2,0.03];
Ymin =[0.1 0.05 0.001];
%-------------------------------------------------------
sname = sprintf('%s/Binned_Exp%d.mat',Prepropath,WHICH_EXP);
X = load(sname,'Fit_Qual','Vector','VectorShuf','ARG');
ARG = X.ARG;

pcrit005 = PCRIT005(METHOD);
pcrit001 = PCRIT001(METHOD);
for par=1:3 % parameter
  for E=1:3
    axes( hAxes{E,par} ); hold on
    % compute group- averages for actual and surrogate data
    Vector = mean(X.Vector{par,E},1);
    VectorShuf = squeeze(mean(X.VectorShuf{par,E},1));
    % consider only frequency dependent models
    Vector = Vector([2:end]);
    VectorShuf = VectorShuf(:,[2:end]);

    % find the correct percentile:
    pc = 100 - 100*(10^(-pcrit005));
    VectorShuf = prctile(VectorShuf,pc);

    plot(ARG.flist,Vector,'-','color',Colorvector(E,:),'LineWidth',1.2);
    plot(ARG.flist,VectorShuf,'-','color',[0.6 0.6 0.6],'LineWidth',1.2);
    issig = find(Vector>VectorShuf);
    plot(ARG.flist(issig),Vector(issig),'.','color',Colorvector(E,:),'MarkerSize',18);

    ylim([Ymin(par) Ymax(par)]);  set(gca,'YTick',[Ymin(par):(Ymax(par)-Ymin(par))/2:Ymax(par)]);  % set(gca,'YTickLabel',[]);
    xlim([1 8]);   set(gca,'XTick',[1:1:8]);    set(gca,'XTickLabel',[]);
    if par==1
      ylabel('Rhythmic component')
      text(2,Ymax(par)*0.9,earname{E},'color',Colorvector(E,:),'FontWeight','bold','FontSize',11)
    end
    if E==3
      xlabel('Frequency [Hz]')
      set(gca,'XTickLabel',[1:1:8]);
    end
    if E==1
      title(sprintf('%s',valnames{par}))
    end
  end % E
end % par



ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',11);



snamef = sprintf('%s/Figure5_effectbinned.jpg',Figurepath);
print('-djpeg','-r600',snamef);


