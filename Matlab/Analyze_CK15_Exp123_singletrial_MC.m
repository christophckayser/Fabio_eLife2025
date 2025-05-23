
clear;
% ---------------------------------------------------------------------
% Monte Carlo simulations of participant samples
% --------------------------------------------------------------------------

ARG =[];
ARG.Modeln = 5;

ARG.flist = [1.2:0.1:4, 4.2:0.2:12]; % frequencies to be tested
ARG.NPerm = 10000; % number of randomizations
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging acvros trials
ARG.Freq_U = 0.5; % frequency of U-term
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset

ARG.Nsim = 5000; % number of draws of participant samples


for WHICH_EXP = 1:4
  WHICH_EXP
   % --------------------------------------------------------------------------
   
  % load data
  defdirsCK15;
  sname = sprintf('%s/ProcessedData_Exp%d.mat',Prepropath,WHICH_EXP);
  x = load(sname);
  Data_all = x.Data_all;
  ARG.VarList = x.ARG.VarList;
  % indictes of reletant variables
  Ind_Freq =  find(strcmp(ARG.VarList,'Freq'));
  Ind_Ear =  find(strcmp(ARG.VarList,'Ear'));
  Ind_SOA =  find(strcmp(ARG.VarList,'SOA'));
  Ind_SOABin =  find(strcmp(ARG.VarList,'SOABin'));
  Ind_Resp = find(strcmp(ARG.VarList,'Resp'));
  Ind_Thr = find(strcmp(ARG.VarList,'SNR'));
  Ind_RT =  find(strcmp(ARG.VarList,ARG.WHICHRT));


  % --------------------------------------------------------------------------
  Nsub = length(Data_all);
  for par=1:3 % parameter
    for E=1:3 % each ear, or combined
      VectorMC{par,E} = zeros(ARG.Nsim,length(ARG.flist)+1);
    end
  end
  for par=1:2
    for E=1:3
      lVector={};
      lARG=[];        lARG = ARG;
      lARG.fast = 1;
      lARG.Modeltype = par;
      parfor R=1:ARG.Nsim
        lVector{R} = zeros(Nsub,length(ARG.flist)+1);
        SubList = ceil(rand(1,Nsub)*Nsub);
        for S=1:Nsub
          Suse = SubList(S);
          tmp = Data_all{Suse}(:,[Ind_RT,Ind_Resp,Ind_SOA,Ind_Thr,Ind_Ear]);
          tmp(:,1) = lARG.rt_transform(tmp(:,1));
          if E<3
            j = find( (tmp(:,5)==E));
          else
            j = [1:size(tmp,1)];
          end
          [~,vs] = local_fitmodelsTrial(tmp(j,[1:3]),lARG);
          lVector{R}(S,:) = squeeze(sqrt(sum(vs([lARG.Modeln-1,lARG.Modeln],:).^2,1)));
        end % S
      end
      
      VectorMC{par,E} = sq(mean(ck_celltomat(lVector),2));
    end
  end

  
  sname = sprintf('%s/SingletrialsMC_Exp%d.mat',Prepropath,WHICH_EXP);
  save(sname,'ARG','VectorMC')
end

return

