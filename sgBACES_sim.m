function [E,D,X] = sgBACES_sim(Y,param_sgBACES) 
  
    zetac = param_sgBACES.zetac;
    zetas = param_sgBACES.zetas;
    nIter = param_sgBACES.nIter;
    Kc = param_sgBACES.mutsha_atoms;
    Ks = param_sgBACES.subspe_atoms;
    bKc = param_sgBACES.mutsha_bases;
    bKs = param_sgBACES.subspe_bases;
    nV = param_sgBACES.nV;
    iM = param_sgBACES.iM; 
    M = param_sgBACES.M;
    lam = param_sgBACES.lambda;
    alpha = param_sgBACES.alpha;
    mu = param_sgBACES.mu;
    mex = param_sgBACES.mex;
 
    N = size(Y{1},1); 
    Dp = dctbases(N,N);
    Dpc = Dp(:,1:bKc);
    Dps = Dp(:,1:bKs);
    iDpc= pinv(Dpc);
    iDps= pinv(Dps);    
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
            R{j} = Y{j}-Dc*Xc;
            if mex == 0
                [As{j},Xs{j}] = my_sBACES_sim(R{j},Dps,iDps,Dpps,As{j},Xs{j},zetas,lam,Ks,nV,alpha,mu);
            else
                [As{j},Xs{j}] = my_sBACES_sim_mex(R{j},Dps,iDps,Dpps,As{j},Xs{j},zetas,lam,Ks,nV,alpha,mu);
            end
            Ds{j} = Dps*As{j};
            E(iter,j+1)  = sqrt(trace((Ds{j}-Dsp{j})'*(Ds{j}-Dsp{j})))/sqrt(trace(Dsp{j}'*Dsp{j}));
        end

        %% Common Level
        for j =1:M
            tmpEc(:,:,j) = Y{j}-Ds{j}*Xs{j};
        end
        Ec = sum(tmpEc,3)./M;
        if mex == 0 
            [Ac,Xc]= my_sBACES_sim(Ec,Dpc,iDpc,Dppc,Ac,Xc,zetac,lam,Kc,nV,alpha,mu);
        else
            [Ac,Xc]= my_sBACES_sim_mex(Ec,Dpc,iDpc,Dppc,Ac,Xc,zetac,lam,Kc,nV,alpha,mu);
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


function [A,X]= my_sBACES_sim(Y,Dp,iDp,Dpp,A,X,zeta,lam,Kc,nV,alpha,mu) 
    D = Dp*A;  
    for j =1:Kc
        X(j,:) = 0;
        F = filter([0 1],1,D);
        Z = (1/(mu+2))*pinv(F)*D*X;
        E = Y-D*X+F*Z;
        xk = D(:,j)'*E;
        xo = spatial(abs(xk),nV,nV);
        thr = zeros(1,size(xk,2));
        for k=1:size(xo,2)
            if abs(xo(k))>alpha*abs(xk(k))
                thr(k) = zeta./abs(xo(k));
            else 
                thr(k) = (zeta./abs(xk(k)))*(1+abs(xo(k)));
            end
        end        
        X(j,:) = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));
        rInd = find(X(j,:));
        if (length(rInd)<1)
            tmp1 = randn(size(Y,1),1);             
            [~,bb]= sort(abs(Dp'*tmp1),'descend');
            ind2 = bb(1:lam);
            A(ind2,j)= (Dp(:,ind2)'*Dp(:,ind2))\Dp(:,ind2)'*tmp1; 
            A(:,j)= A(:,j)./norm(Dp*A(:,j));
            D(:,j) = Dp*A(:,j);
            tmp2 = D(:,j)'*Y;
            thr = zeta./abs(tmp2);
            X(j,:) = sign(tmp2).*max(0, bsxfun(@minus,abs(tmp2),thr/2));
            continue
        end
        mag = X(j,rInd)*X(j,rInd)';
        tmp3 = D(:,j)+(1/mag)*E(:,rInd)*X(j,rInd)';
        [~,bb]= sort(abs(Dp'*tmp3),'descend');
        ind = bb(1:lam);
        A(ind,j)= (Dp(:,ind)'*Dp(:,ind))\Dp(:,ind)'*tmp3;
        A(:,j) = A(:,j)./norm(Dp*A(:,j));
        D(:,j) = Dp*A(:,j);
    end
    A = Pruning(Dp,A,X,Y,iDp,Dpp);
end


function A = Pruning(Dp,A,X,Y,invDp,Dpp)
    % idea adopted from M. Aharon et.al, "K-SVD: An algorithm for designing overcomplete dictionaries for SR," in IEEE Transactions on Signal Processing,
    % vol. 54, no. 11, pp. 4311-4322, Nov. 2006, doi: 10.1109/TSP.2006.881199.
    TCs_thresh = 0.99;
    SMs_thresh = 3;
    E = sum((Y-Dp*A*X).^2,1);
    R = A'*Dpp*A;
    R = R-diag(diag(R));
    for jj=1:size(A,2)
        if max(R(jj,:))>TCs_thresh  || length(find(abs(X(jj,:))>1e-5))<=SMs_thresh
            [~,ind]=max(E);
            E(ind(1))=0;
            A(:,jj) = invDp*Y(:,ind(1));
            A(:,jj)= A(:,jj)./norm(Dp*A(:,jj));
            R = A'*Dpp*A;
            R = R-diag(diag(R));
        end
    end
end

function xo = spatial(x,nV1,nV2)
    x = reshape(x,nV1,nV2);
    xr = padarray(x,[1 1],0,'both');
    TT = conv2(xr,[1 1 1;1 1 1;1 1 1],'valid');
    xo = reshape(TT,1,nV1*nV2);
end















