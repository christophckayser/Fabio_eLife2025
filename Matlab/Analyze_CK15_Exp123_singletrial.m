
clear;
% ---------------------------------------------------------------------
% Analysis of behavior time course for CK 15 Exp 1-4
% Single trial data against SOA using linear models
% --------------------------------------------------------------------------
ARG =[];
ARG.Modeln = 5;

ARG.flist = [1.2:0.1:4, 4.2:0.2:12]; % frequencies to be tested
ARG.NPerm = 10000; % number of randomizations
ARG.rt_transform = @(x) sqrt(x); % transformation of RTs before averaging acvros trials
ARG.Freq_U = 0.5; % frequency of U-term
ARG.WHICHRT = 'RT_Target'; % RT_Target RT_Onset



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
      Fit_Qual{par,E} = zeros(Nsub,length(ARG.flist)+1);
      Vector{par,E} = zeros(Nsub,length(ARG.flist)+1);
      VectorShuf{par,E} = zeros(Nsub,ARG.NPerm,length(ARG.flist)+1);
    end
  end
  for par=1:2
    for E=1:3
      lVector={};
      lVectorShuf={};
      lFit_Qual={};
      parfor S=1:Nsub
        lARG=[];        lARG = ARG;
        lFit_Qual{S} = zeros(1,length(lARG.flist)+1);
        lVector{S} = zeros(1,length(ARG.flist)+1);
        lVectorShuf{S} = zeros(lARG.NPerm,length(ARG.flist)+1);

        % --------------------------------------------------------------------------
        % actual data
        lARG.fast = 0;
        lARG.Modeltype = par;
        tmp = Data_all{S}(:,[Ind_RT,Ind_Resp,Ind_SOA,Ind_Thr,Ind_Ear]);
        tmp(:,1) = lARG.rt_transform(tmp(:,1));
        if E<3
          j = find( (tmp(:,5)==E));
        else
          j = [1:size(tmp,1)];
        end
        [lFit_Qual{S},vs] = local_fitmodelsTrial(tmp(j,[1:3]),lARG);
        lVector{S}(1,:) = squeeze(sqrt(sum(vs([lARG.Modeln-1,lARG.Modeln],:).^2,1)));

        % --------------------------------------------------------------------------
        %  same for surrogate data
        lARG.fast = 1; % don't compute LL
        for rep=1:lARG.NPerm
          R = randperm(size(tmp,1));
          tmp(:,[1,2]) = tmp(R,[1,2]);
          if E<3
            j = find( (tmp(:,5)==E));
          else
            j = [1:size(tmp,1)];
          end
          [~,vs] = local_fitmodelsTrial(tmp(j,[1:3]),lARG);
          lVectorShuf{S}(rep,:) =  single(squeeze(sqrt(sum(vs([lARG.Modeln-1,lARG.Modeln],:).^2,1))));
        end
      end % parfor S

      Vector{par,E} = ck_celltomat(lVector);
      VectorShuf{par,E} = ck_celltomat(lVectorShuf);
      Fit_Qual{par,E} = ck_celltomat(lFit_Qual);
    end
  end

  
  sname = sprintf('%s/Singletrials_Exp%d.mat',Prepropath,WHICH_EXP);
  save(sname,'Fit_Qual','ARG','Vector','VectorShuf')
end

return

