function [E,D,X] = ShSSDL_ADMM_sim(Y,param_ShSSDL)
    zetac = param_ShSSDL.zetac;
    zetas = param_ShSSDL.zetas; 
    nIter1 = param_ShSSDL.nIter_outer;
    nIter2 = param_ShSSDL.nIter_inner;
    alpha = param_ShSSDL.coherence;
    Kc = param_ShSSDL.mutsha_atoms;
    Ks = param_ShSSDL.subspe_atoms;
    spams = param_ShSSDL.spams;
    paramc.L = zetac; 
    paramc.eps=0.001; 
    params.L = zetas; 
    params.eps=0.001; 
    iM = param_ShSSDL.iM;
    M = param_ShSSDL.M;
    N = size(Y{1},1); 

    for j=1:M
        if j == iM
            Dc = randn(N,Kc);% Y{j}(:,1:Kc);  
            Dc = Dc*diag(1./sqrt(sum(Dc.*Dc)));
            Xc = zeros(Kc,size(Y{j},2)); 
        end
        Ds{j} = randn(N,Ks); %Y{j}(:,1:Ks);  
        Ds{j} = Ds{j}*diag(1./sqrt(sum(Ds{j}.*Ds{j})));
        Xs{j} = zeros(Ks,size(Y{j},2));
    end   

    for iter = 1:nIter1
        Dcp = Dc;
        for j =1:M
            Dsp{j} = Ds{j};
        end       

        %% sparse coding 
        for iter2 = 1:nIter2
            %%% Subjectwise
            for j=1:M
                Eg{j} = Y{j}-Dc*Xc;
                Xs{j} = full(mexOMP(Eg{j},Ds{j},params));
            end

            %%% Common       
            for j=1:M
                tmpEc(:,:,j) = Y{j}-Ds{j}*Xs{j};
            end
            Ec = sum(tmpEc,3)./M;
            Xc = full(mexOMP(Ec,Dc,paramc));
        end

        %% Dictionary update

        St.nSub=M; St.eps = 10^-4; St.AlgoType=2; St.Kmax = 40; %st.kmax is number of iterations
        St.eta = alpha; %0.5; %coherence parameter
        St.scale_mu = 2;  St.ini_mu = 10^-4;  St.max_mu = 10^10; %ADMM parameters 
        [Dc,Ds] = DictUpdate2(Y,Dc,Ds,Xc,Xs,St);
        E(iter,1)= sqrt(trace((Dc-Dcp)'*(Dc-Dcp)))/sqrt(trace(Dcp'*Dcp));
        for j=1:M
            E(iter,j+1) = sqrt(trace((Ds{j}-Dsp{j})'*(Ds{j}-Dsp{j})))/sqrt(trace(Dsp{j}'*Dsp{j}));
        end


    end
    
    Dss=[];
    Xss=[];
    for i =1:M
        Dss = [Dss Ds{i}];
        Xss = [Xss;Xs{i}];
    end

    D = [Dc Dss];
    X = [Xc;Xss];
end



%% Dictionary Update Step Using ADMM
function [Dict_0,Dict] = DictUpdate2(Data,Dict_0,Dict,X_0,X,St)
    M = St.nSub;    max_mu = St.max_mu;     scale_mu = St.scale_mu;   ini_mu = St.ini_mu;
    Z0 = zeros(size(Dict_0));   mu = ini_mu;    
    Zi = zeros(size(Dict{1}));  W0 = Z0;    Wi = Zi; 
    if St.AlgoType == 1
        DictCommon;
        DictIndividual;
    elseif St.AlgoType == 2
        DictIndividual;
        DictCommon;
    end

    % Inner Functions
    function DictCommon
        T = 0;  A0 = [];
        Z0 = zeros(size(Dict_0));   mu = ini_mu; W0 = Z0;
        for i = 1:M    
            A0 = [A0,Dict{i}];
            T = T + (Data{i} - Dict{i} * X{i});
        end
        DD = Dict_0;    T = T./M;   
        A00 = A0*A0';   %rho = 1/norm(full(X_0)*full(X_0)');
        for k = 1:St.Kmax
%             while norm((DD - Z0),'fro') > 0.01
                DD = (mu*Z0 + T*X_0' - W0)/(X_0*X_0' + mu*eye(size(DD,2)));  DD = normc(DD);
                Z0 = (mu*eye(size(DD,1))+ 2*St.eta*A00)\(W0 + mu*DD); Z0 = normc(Z0);
%             end
            W0 = W0 + mu*(DD-Z0);
            mu = min(scale_mu*mu,max_mu);
%             display(norm((DD - Z0),'fro'))
            if norm((DD - Z0),'fro') < St.eps;  break;  end
        end
        Dict_0 = DD;
    end

    function DictIndividual
        for i = 1:M        
            Zi = zeros(size(Dict{i}));   mu = ini_mu; Wi = Zi;
            Ai = Dict_0;
            for h = 1:M   
                if h == i; continue; end 
                Ai = [Ai,Dict{h}];   
            end;
            Gi = Data{i} - Dict_0*X_0;
            DD = Dict{i};   XX = X{i};
            Aii = Ai*Ai';   % rho = 1/norm(full(XX)*full(XX)');
            for k = 1:St.Kmax
%                 while norm((DD - Zi),'fro') > 0.01
                    DD = (mu*Zi + Gi*XX' - Wi)/(XX*XX' + mu*eye(size(DD,2)));  DD = normc(DD);
                    Zi = (mu*eye(size(DD,1))+ 2*St.eta*Aii)\(Wi + mu*DD); Zi = normc(Zi);
%                 end
                Wi = Wi + mu*(DD-Zi);
                mu = min(scale_mu*mu,max_mu);
            
                if norm((DD - Zi),'fro') < St.eps;    break;         end
            end
            Dict{i} = (DD);
        end        
    end
end








