clear 
close all
clc

load Simtb_data

%% general parameters
stdv(1,:) = [0.0 0.0 0.0 0.3 0.3 0.3 0.6 0.6 0.6 0.9 0.9 0.9]; 
stdv(2,:) = repmat([0.06 0.12 0.18],1,4); %0.06 0.12 0.18
M = 8;
N = 230;
nV = 3600;
nCC = 4; % number of common components
nIter = 20;
Kc = 12;
Ks = 4;  
thresh = 0;

%% Selecting the subject whose TC and SMs to be used as the ground truth
iM = 1; % papers results are according to subject = 1
gSM = reshape(allSM(iM,1: nCC,:),nCC,nV);
gTC = reshape(zscore(allTC(iM,:,1:nCC)),N,nCC);

SM = cell(1,M);
TC = cell(1,M);
for sub=1:M
    SM{sub} = reshape(allSM(sub,[1:nCC nCC+sub],:),nCC+1,nV);
    TC{sub} = reshape(zscore(allTC(sub,:,[1:nCC nCC+sub])),N,nCC+1);
    gSM = [gSM; SM{sub}(nCC+1,:)];
    gTC = [gTC TC{sub}(:,nCC+1)];
end

%% algorithm specific parameters
param_CODL.spa = 1.5; 
param_CODL.iter = 3*nIter;
param_CODL.batsize = nV;
param_CODL.atoms = Kc+Ks;
param_CODL.tol = 1e-9;
param_CODL.spams = 0; % set it to 1 if you want to use SPAMS toolbox dictionary  

param_ShSSDL.zetac = 2; 
param_ShSSDL.zetas = 1; 
param_ShSSDL.nIter_outer = nIter;
param_ShSSDL.nIter_inner = 3;
param_ShSSDL.coherence = 2; 
param_ShSSDL.mutsha_atoms = Kc;
param_ShSSDL.subspe_atoms = Ks;
param_ShSSDL.cd_thresh = thresh;
param_ShSSDL.spams = 1; % set it to 0 if you do not have SPAMS toolbox installed, the code will instead use matlab implementation of OMP, which is very slow
param_ShSSDL.iM = iM;
param_ShSSDL.M = M;

param_rgBACESW.zetac = 2; %4
param_rgBACESW.zetas = 14; %12
param_rgBACESW.nIter = nIter;
param_rgBACESW.mutsha_atoms = Kc;
param_rgBACESW.subspe_atoms = Ks;
param_rgBACESW.mutsha_bases = 120;
param_rgBACESW.subspe_bases = 120;
param_rgBACESW.Atol = 1e-5;
param_rgBACESW.iM = iM;
param_rgBACESW.M = M;
param_rgBACESW.mex = 0;

param_sgBACES.zetac = 11;
param_sgBACES.zetas = 125; 
param_sgBACES.nIter = nIter;
param_sgBACES.mutsha_atoms = Kc;
param_sgBACES.subspe_atoms = Ks;
param_sgBACES.mutsha_bases = 120;
param_sgBACES.subspe_bases = 120;
param_sgBACES.nV = sqrt(nV);
param_sgBACES.iM = iM;
param_sgBACES.M = M;
param_sgBACES.lambda = 105;
param_sgBACES.alpha = 1.5; 
param_sgBACES.mu = 3; 
param_sgBACES.mex = 0; % 1 for running the mex functions instead of matlab functions 
                      % paper's results are according to mex = 1; 

%% seed
% rng('default')
% rng(50,'philox') 
Y = cell(1,M);
   
%% Simulation
nAlg = 5;
Trials = 5;
RC = 2; % recovery cases
S = size(stdv,2);
R = [1 1 0 0 0]; % vector needed for correlation calculations: subject-specific or group
tStart_main = tic;
for t = 1:Trials
    fprintf(2,'Trial #%.i\n',t)
    for j =1:S
        for sub=1:M
            Y{sub} = (TC{sub}+stdv(1,j)*randn(N,nCC+1))*(SM{sub}+stdv (2,j)*randn(nCC+1,nV)); 
            Y{sub} = Y{sub}-repmat(mean(Y{sub}),size(Y{sub},1),1); 
        end
    
        % sICA
        tStart = tic;
        [rTC{1},rSM{1}] = sICA_sim(Y,M,Kc+Ks,Kc+Ks); 
        tEnd(t,j,1) = toc(tStart);
        fprintf('spatial group ICA took %f seconds\n',  rem(tEnd(t,j,1),60));

        % CODL
        tStart = tic; 
        [E1,rTC{2},rSM{2}] = CODL_sim(Y,M,param_CODL); 
        tEnd(t,j,2) = toc(tStart);
        fprintf('CODL took %f seconds\n',  rem(tEnd(t,j,2),60));       

        % ShSSDL
        tStart = tic; 
        [E2,rTC{3},rSM{3}] = ShSSDL_ADMM_sim(Y,param_ShSSDL);
        tEnd(t,j,3) = toc(tStart);
        fprintf('ShSSDL took %f seconds\n', rem(tEnd(t,j,3),60));       
    
        % rgBACESW   
        tStart = tic;     
        [~,rTC{4},rSM{4}]= rgBACESW_sim(Y,param_rgBACESW);  
        tEnd(t,j,4) = toc(tStart);
        fprintf('sgBACES took %f seconds\n',  rem(tEnd(t,j,4),60));       

        % sgBACES   
        tStart = tic;     
        [~,rTC{5},rSM{5}]= sgBACES_sim(Y,param_sgBACES);  
        tEnd(t,j,5) = toc(tStart);
        fprintf('sgBACES took %f seconds\n',  rem(tEnd(t,j,5),60));       
        
        for k=1:nAlg
            for l=1:RC
                if l==1
                    [srTC{k},srSM{k}] = sort_TSandSM_temporal(gTC,gSM,rTC{k},rSM{k},Kc,M,Ks,nCC,R(k));
                else
                    [srTC{k},srSM{k}] = sort_TSandSM_spatial(gTC,gSM,rTC{k},rSM{k},Kc,M,Ks,nCC,R(k));
                end
                TCcorr(:,k,j,t) = abs(diag(corr(gTC,srTC{k}))); 
                mCorr(l,1,k,j,t) = mean(TCcorr(:,k,j,t));
                SMcorr(:,k,j,t) = abs(diag(corr(abs(gSM'),abs(srSM{k}')))); 
                mCorr(l,2,k,j,t) = mean(SMcorr(:,k,j,t));
            end
        end
        
        fprintf('sICA mCTC_t = %.4f    ', mCorr(1,1,1,j,t)); 
        fprintf('CODL mCTC_t = %.4f    ', mCorr(1,1,2,j,t)); 
        fprintf('ShSSDL mCTC_t = %.4f    ', mCorr(1,1,3,j,t));        
        fprintf('rgBACESW mCTC_t = %.4f   ', mCorr(1,1,4,j,t)); 
        fprintf('sgBACES mCTC_t = %.4f \n', mCorr(1,1,5,j,t));
        
        fprintf('sICA mCSM_t = %.4f    ', mCorr(1,2,1,j,t)); 
        fprintf('CODL mCSM_t = %.4f    ', mCorr(1,2,2,j,t)); 
        fprintf('ShSSDL mCSM_t = %.4f    ', mCorr(1,2,3,j,t));      
        fprintf('rgBACESW mCSM_t = %.4f   ', mCorr(1,2,4,j,t)); 
        fprintf('sgBACES mCSM_t = %.4f \n', mCorr(1,2,5,j,t));   

        fprintf('sICA mCTC_s = %.4f    ', mCorr(2,1,1,j,t)); 
        fprintf('CODL mCTC_s = %.4f    ', mCorr(2,1,2,j,t)); 
        fprintf('ShSSDL mCTC_s = %.4f    ', mCorr(2,1,3,j,t));      
        fprintf('rgBACESW mCTC_s = %.4f   ', mCorr(2,1,4,j,t)); 
        fprintf('sgBACES mCTC_s = %.4f \n', mCorr(2,1,5,j,t));

        fprintf('sICA mCSM_s = %.4f    ', mCorr(2,2,1,j,t)); 
        fprintf('CODL mCSM_s = %.4f    ', mCorr(2,2,2,j,t)); 
        fprintf('ShSSDL mCSM_s = %.4f    ', mCorr(2,2,3,j,t));      
        fprintf('rgBACESW mCSM_s = %.4f   ', mCorr(2,2,4,j,t)); 
        fprintf('sgBACES mCSM_s = %.4f \n', mCorr(2,2,5,j,t));
        fprintf('\n') 

    end
end
mean(mean(tEnd(:,:,:)))
tEnd_main = toc(tStart_main);
fprintf('%d minutes and %f seconds\n', floor(tEnd_main/60), rem(tEnd_main,60));


%%
figure;
NV= {'0.0, 0.06' '0.0, 0.12' '0.0, 0.18' '0.3, 0.06' '0.3, 0.12' '0.3, 0.18' '0.6, 0.06' '0.6, 0.12' '0.6, 0.18' '0.9, 0.06' '0.9, 0.12' '0.9, 0.18'}; 
TT = mean(mCorr,5);
Markers = {'+','o','*','x','v','d','^','s','>','<','p','h'};
A = {'Source retrieval with respect to GT-TCs','Source retrieval with respect to GT-TCs',...
     'Source retrieval with respect to GT-SMs','Source retrieval with respect to GT-SMs'};
B = {'Mean of mcTC_t','Mean of mcSM_t','Mean of mcTC_s','Mean of mcSM_s'};
c = lines(nAlg);
count = 0;
for k= 1:RC
    for i =1:2
        for j =1:nAlg
            subplot(2,2,2*(k-1)+i);
            tmp1(:,j) = reshape(TT(k,i,j,:),S,1);
            plot(tmp1(:,j), 'Marker',Markers{j},'Color', c(j,:),"MarkerSize",5,'MarkerIndices',1:length(tmp1(:,j))); hold on;
        end
        grid on; ax = gca; ax.XTick = 1:12; ax.XTickLabel = ''; set(gca,'Ytick',min(min(tmp1))-0.025:0.05:1)
        ylabel(B{2*(k-1)+i}, 'FontSize', 12)
        xlabel('Spatiotemporal Noise Std')
        title (A{2*(k-1)+i})
        set(gca,'xtick',1:12,'xticklabel',NV)
    end
end
legend('sICA','CODL','ShSSDL','rgBACESW','sgBACES')