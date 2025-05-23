
clear; close all;

% Figure 9 showing with eye data for experiment 4. only for those trials with RT > 0.9 sec
% to avoid confounds from motor action

clear;

defdirsCK15

sname = sprintf('%s/Data_for_eyefigure.mat',Prepropath);
load(sname)

tstr = {'X-pos mean','X-pos SD','Y-pos mean','Y-pos SD','Pupil diameter','# Saccades'};

figure('Position',[94         233        1190         653])
nRows = 2 ; nCols = 3 ;
% - Create grid of axes.
[blx, bly] = meshgrid( 0.1:0.9/nCols:0.82, 0.1:0.85/nRows:0.8 ) ;
bly = bly([nRows:-1:1],:);
blx(:,[2]) = blx(:,[2])-0.022;

hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, 0.68/nCols, 0.61/nRows] ), blx, bly, 'UniformOutput', false ) ;


ranges(1,:) = [-0.5 0.5];
ranges(2,:) = [0 5];
ranges(3,:) = [-0.5 0.5];
ranges(4,:) = [0 5];
ranges(5,:) = [-0.2 0.2];
ranges(6,:) = [0 0.08];
vpos = [1,2,3,4,5,6];

ylab = {'deg','deg','deg','deg','Z','per Trial'}
for v=1:6
  axes( hAxes{vpos(v)} ); hold on
  tmp = sq(AllData(:,v,:));

  % time 50 is 0ms
  if sum(v==[1,3,5])
    tmp = tmp-tmp(:,50);
  end
  plot(taxf,tmp,'color',ones(1,3)*0.4);
  hold on
  plot(taxf,mean(tmp),'color',[0 0 0 ],'LineWidth',2.5)
  xlim([-0.4 1]); set(gca,'XTick',[-0.2:0.2:1],'XTickLabel',[]);
  title(tstr{v})
  ylim([ranges(v,1) ranges(v,2)])
  if rem(v,2)==0
    xlabel('Time [s]');
    set(gca,'XTick',[-0.2:0.2:1],'XTickLabel',[-0.2:0.2:1]);
  end
  ylabel(ylab{v})
end


% -----------------------------------------
% compute mean sacade date in tiem range of 0 to 1.2s
int = find( (taxf>=0).*(taxf<=1.2));
Saccrate = sum(AllData(:,6,int),3)
[mean(Saccrate),sem(Saccrate)]


% -----------------------------------------
% add analysis labels
axes( hAxes{1} );
text(-0.67,0.65,'A','FontSize',20,'FontWeight','bold');

axes( hAxes{2} );
text(-0.67,5.8,'B','FontSize',20,'FontWeight','bold');

axes( hAxes{5} );
text(-0.67,0.26,'C','FontSize',20,'FontWeight','bold');

axes( hAxes{6} );
text(-0.67,0.09,'D','FontSize',20,'FontWeight','bold');

ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',11);


snamef = sprintf('%s/Figure9_eyes.jpg',Figurepath);
print('-djpeg','-r600',snamef);

