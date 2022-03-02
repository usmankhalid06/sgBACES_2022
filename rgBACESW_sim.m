function [E,D,X] = rgBACESW_sim(Y,param_rgBACESW) 
  
    zetac = param_rgBACESW.zetac;
    zetas = param_rgBACESW.zetas;
    nIter = param_rgBACESW.nIter;
    Kc = param_rgBACESW.mutsha_atoms;
    Ks = param_rgBACESW.subspe_atoms;
    bKc = param_rgBACESW.mutsha_bases;
    bKs = param_rgBACESW.subspe_bases;
    tol = param_rgBACESW.Atol;
    iM = param_rgBACESW.iM; 
    M = param_rgBACESW.M;
    mex = param_rgBACESW.mex;
    
    N = size(Y{1},1); 
    Dp = dctbases(N,N);
    Dpc = Dp(:,1:bKc);
    Dps = Dp(:,1:bKs);
    invDpc= pinv(Dpc);
    invDps= pinv(Dps);
    Dppc= Dpc'*Dpc;
    Dpps= Dps'*Dps;

    for j=1:M
        if j ==iM
            Dc = Dp(:,1:Kc);
            Ac = eye(bKc,Kc); 
            Xc = zeros(Kc,size(Y{j},2)); 
        end 
        Ds{j} = Dp(:,1:Ks);
        As{j} = eye(bKs,Ks);
        Xs{j} = zeros(Ks,size(Y{j},2));
    end   
    

    for iter = 1:nIter  
        Dcl = Dc;
        for j =1:M
            Dsp{j} = Ds{j};
        end

        %% Subject Level
        for j=1:M
            ADs2 = [];
            for i =[1:j-1 j+1:M]
                ADs2 = [ADs2 Dps*As{i}];
            end
            ADs2 = [Dpc*Ac ADs2];
            R{j} = Y{j}-Dc*Xc;
            if mex == 0
                [As{j},Xs{j}] = my_rBACESW(R{j},Dps,invDps,Dpps,As{j},Xs{j},zetas,tol,Ks,bKs);
            else
                [As{j},Xs{j}] = my_rBACESW_mex(R{j},Dps,invDps,Dpps,As{j},Xs{j},zetas,tol,Ks,bKs);
            end
            Ds{j} = Dps*As{j};
            E(iter,j+1)  = sqrt(trace((Ds{j}-Dsp{j})'*(Ds{j}-Dsp{j})))/sqrt(trace(Dsp{j}'*Dsp{j}));
        end

        %% Common Level
        ADs1 = [];
        for j =1:M
            ADs1 = [ADs1 Dps*As{j}];
            tmpEc(:,:,j) = Y{j}-Ds{j}*Xs{j};
        end
        Ec = sum(tmpEc,3)./M;
        if mex == 0
            [Ac,Xc]= my_rBACESW(Ec,Dpc,invDpc,Dppc,Ac,Xc,zetac,tol,Kc,bKc);
        else
            [Ac,Xc]= my_rBACESW_mex(Ec,Dpc,invDpc,Dppc,Ac,Xc,zetac,tol,Kc,bKc);
        end
        Dc = Dpc*Ac;
        E(iter,1) = sqrt(trace((Dc-Dcl)'*(Dc-Dcl)))/sqrt(trace(Dcl'*Dcl));

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
   
function D = dctbases(N,K)
    D = zeros(N,K);
    for k = 1:K
        d = cos(pi*(0:N-1)*(k-1)/K)';
        if k>1
            d = d-mean(d);
        else
            d = d-0;
        end
        D(:,k) = d/norm(d);
    end
end



function [A,X]= my_rBACESW(Y,Dp,invDp,Dpp,A,X,zeta,tol,Kc,bKc)
    D = Dp*A; 
    for j =1:Kc
        Ao= zeros(bKc,1);
        X(j,:) = 0;
        E = Y-D*X;
        xk = D(:,j)'*E;
        thr = zeta./abs(xk);  
        while norm(A(:,j)-Ao)/norm(Ao)>tol
            Ao = A(:,j);
            tmp = (Dp*A(:,j))'*E;
            X(j,:) = sign(tmp).*max(0, bsxfun(@minus,abs(tmp),thr/2));
            mag = X(j,:)*X(j,:)';
            A(:,j)= invDp*(D(:,j)+(1/mag)*(E*X(j,:)'));
            A(:,j) = A(:,j)./norm(Dp*A(:,j));
        end
        D(:,j) = Dp*A(:,j);
        rInd = find(X(j,:));
        if (length(rInd)<1)
            tmp1 = randn(size(Y,1),1);   
            A(:,j) = invDp*tmp1;
            A(:,j)= A(:,j)./norm(Dp*A(:,j));
            D(:,j) = Dp*A(:,j);
            tmp = (D(:,j))'*Y;
            thr = zeta./abs(tmp);            
            X(j,:) = sign(tmp).*max(0, bsxfun(@minus,abs(tmp),thr/2));
        end
    end
    A = Pruning(Dp,A,X,Y,invDp,Dpp);
end


function A = Pruning(Dp,A,X,Y,invDp,Dpp)
    % idea adopted from M. Aharon et.al, "K-SVD: An algorithm for designing overcomplete dictionaries for SR," in IEEE Transactions on Signal Processing, 
    % vol. 54, no. 11, pp. 4311-4322, Nov. 2006, doi: 10.1109/TSP.2006.881199.
    TCs_thresh = 0.99;
    SMs_thresh = 3;
    Er = sum((Y-Dp*A*X).^2,1); 
    R = A'*Dpp*A; 
    R = R-diag(diag(R));
    for jj=1:size(A,2)
        if max(R(jj,:))>TCs_thresh  || length(find(abs(X(jj,:))>1e-5))<=SMs_thresh 
            [~,ind]=max(Er);
            Er(ind(1))=0;
            A(:,jj) = invDp*Y(:,ind(1));
            A(:,jj)= A(:,jj)./norm(Dp*A(:,jj));
            R = A'*Dpp*A; 
            R = R-diag(diag(R));
        end
    end
end






