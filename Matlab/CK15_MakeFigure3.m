
clear;

% Visualize the data for experiments 1-3.
% Show d-prime, bias and RT for time-binne
% Figure 3

% arguments
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging over trials
ARG.WHICHRT = 'RT_Target'; % which RTs to show. 

% Binning of SOA in ms
ARG.BINS = [0:60:1200]./1000; %  20 time bins, each 60ms, max frequency is 8.3 Hz 
nbins = length(ARG.BINS)-1;
% colors for display
Colorvector(1,:) = [250,150,40]/250; %
Colorvector(2,:) = [180,40,40]/250; %
Colorvector(3,:) = [40,150,250]/250; %
Colorvector(4,:) = [40,40,180]/250; %

Colorvector2 = [0.3,0.3,0.3];
Colorvector2 = repmat(Colorvector2,[10,1]);

% make sure directories exist.
defdirsCK15;

figure('Position',[    102         173        1489         769])
% ------------------------------------------------------
nRows = 3 ; nCols = 6 ;
% - Create grid of axes.
[blx, bly] = meshgrid( 0.14:0.9/nCols:0.9, 0.1:0.8/nRows:0.8 ) ;
bly = bly([nRows:-1:1],:);
blx(:,[2]) = blx(:,[2])-0.07;
blx(:,[3]) = blx(:,[3])-0.14;

blx(:,[4]) = blx(:,[4])-0.18;
blx(:,[5]) = blx(:,[5])-0.14;
blx(:,[6]) = blx(:,[6])-0.1;

ColSize = [0.05 0.05 0.08 0.16 0.16 0.16 ];
ColSize = repmat(ColSize,[3,1]);
hAxes = arrayfun( @(x,y,z) axes( 'Position', [x, y, z, 0.66/nRows] ), blx, bly, ColSize, 'UniformOutput', false ) ;


% ------------------------------------------------------------------
% loop experiments

for WHICH_EXP = 1:3

  sname = sprintf('%s/ProcessedData_Exp%d.mat',Prepropath,WHICH_EXP);
  x = load(sname);
  ARG.VarList = x.ARG.VarList;

  % ------------------------------------------------------------
  % get data and bin
  for S=1:length(x.Data_all)
    Behav(S) = compute_binned_behavior(x.Data_all{S},ARG);
  end
  % combine data across participants. converts formats
  Behav = local_combinebehavpart(Behav);
  ns = size(Behav.pc,1);

  % ------------------------------------------------------------
  % Display data

  % behavior independent of SOA
  axes( hAxes{WHICH_EXP,1} )

  % d prime
  ckmeanplotcompact(squeeze(Behav.dp),1,1,1,Colorvector);
  set(gca,'XTick',[1:2],'XTickLabel',{});
  axis([0.7 2.6 0 3.5]);
  if WHICH_EXP==1
    title('d-prime')
  elseif WHICH_EXP==3
    set(gca,'XTick',[1:2],'XTickLabel',{'L','R'});
  end
  % bias
  axes( hAxes{WHICH_EXP,2} )
  ckmeanplotcompact(squeeze(Behav.crit),1,1,1,Colorvector);
  set(gca,'XTick',[1:2],'XTickLabel',{});
  axis([0.7 2.6 -1.5 1.5])
  if WHICH_EXP==1
    title('bias')
  elseif WHICH_EXP==3
    set(gca,'XTick',[1:2],'XTickLabel',{'L','R'});
  end

  % rt
  axes( hAxes{WHICH_EXP,3} )
  tmp = reshape(Behav.rt,[ns,4]);
  ckmeanplotcompact(tmp,1,1,1,Colorvector);
  set(gca,'XTick',[1:4],'XTickLabel',{});
  axis([0.7 4.8 0.5 1.5])
  if WHICH_EXP==1
    title('rt ')
  elseif WHICH_EXP==3
    set(gca,'XTick',[1:4],'XTickLabel',{'Lf1','Lf2','Rf1','Rf2'});
  end


  % -------------------------------------------------------------------
  % time course for ears combined
  E=1;

  % dp per time bin
  axes( hAxes{WHICH_EXP,4} )
  tmp = (squeeze(Behav.dpSOA(:,E,:))')';
  ckmeanplotcompact(tmp,1,1,1,Colorvector2);
  axis([0.5 nbins+1 0 3.5]); set(gca,'XTick',[1:5:21])
  if WHICH_EXP==1
    title('d-prime')
  elseif WHICH_EXP==3
    xlabel('time'); set(gca,'XTickLabel',ARG.BINS([1:5:21]));
  end

  % bias per time bin
  axes( hAxes{WHICH_EXP,5} )
  tmp = (squeeze(Behav.critSOA(:,E,:))')';
  ckmeanplotcompact(tmp,1,1,1,Colorvector2);
  axis([0.5 nbins+1 -1.5 1.5]); set(gca,'XTick',[1:5:21],'YTick',[-1.5,0,1.5])
  if WHICH_EXP==1
  title('bias')
  elseif WHICH_EXP==3
    xlabel('time');  set(gca,'XTickLabel',ARG.BINS([1:5:21]));
  end

  % rt per time bin
  axes( hAxes{WHICH_EXP,6} )
  tmp = squeeze(Behav.rtSOA(:,E,:));
  ckmeanplotcompact(tmp,1,1,1,Colorvector2);
  axis([0.5 nbins+1 0.5 1.5]); set(gca,'XTick',[1:5:21],'YTick',[0.5:0.5:1.5])
  if WHICH_EXP==1
    title('rt')
  elseif WHICH_EXP==3
    xlabel('time');  set(gca,'XTickLabel',ARG.BINS([1:5:21]));
  end

end


% figure cosmetics
ckfigure_setall(gcf,'TickLength',[0.02 0.02]);
ckfigure_setall(gcf,'Box','Off');
ckfigure_setall(gcf,'FontSize',12);

% add experimentl labels
for WHICH_EXP=1:3
  axes( hAxes{WHICH_EXP,1} );
  text(-0.5,1.5,sprintf('Exp %d',WHICH_EXP),'FontSize',16,'Rotation',90);

end

axes( hAxes{1,1} );
text(-1,4.2,'A','FontSize',20,'FontWeight','bold');

axes( hAxes{1,4} );
text(-2,4.2,'B','FontSize',20,'FontWeight','bold');




snamef = sprintf('%s/Figure3_maindata.jpg',Figurepath);
print('-djpeg','-r300',snamef);


