function [E,Dg,Xg] = CODL_sim(Y,nS,param_CODL) 

    zeta = param_CODL.spa;
    nIter = param_CODL.iter;
    bs = param_CODL.batsize;
    dic_size = param_CODL.atoms;
    tol = param_CODL.tol;
    spams = param_CODL.spams; 

    FY = [];
    for j=1:nS
        FY= [FY;Y{j}];
    end
    N = size(Y{1},1);

    param.K=dic_size;  % learns a dictionary with 100 elements
    param.lambda=zeta;
    param.numThreads=4; % number of threads
    param.batchsize=bs;
    param.iter=nIter;  % let us see what happens after 1000 iterations.

    if spams ==1
        tmpDg = mexTrainDL(FY,param);
    else
        Di = FY(:,1:dic_size); 
        [tmpDg,~,E] = my_ODL(FY,Di,zeta,tol,nIter);
    end
    
    Xg=mexLasso(FY,tmpDg,param);
    
    for i=1:dic_size
        [Dg(:,i),~,~] = svds(reshape(tmpDg(:,i),N,nS),1);
    end

end



function [D, X, Err] = my_ODL(Y, Di, lambda, tol, nIter)
    D = Di;
    X = zeros(size(D,2), size(Y,2));
	iter = 0;
    param.mode = 2; 
    param.lambda = lambda; 
	while (iter < nIter)
        Dold = D;
		iter = iter + 1;
        X = mexLasso(Y, D, param);
        XX = X*X';
        YX = Y*X';
        Dp = D;
        iter2 = 0;
        while (iter2 < nIter)
            iter2 = iter2 + 1;
            for i = 1: size(D,2)
                if(XX(i,i) ~= 0)
                    a = 1.0/XX(i,i) * (YX(:,i) - D*XX(:, i)) + D(:,i);
                    D(:,i) = a/(max( norm(a,2),1));
                end
            end
            if (norm(D - Dp, 'fro')/numel(D) < tol)
                break;
            end
            Dp = D;
        end
        Err(iter)= sqrt(trace((D-Dold)'*(D-Dold)))/sqrt(trace(Dold'*Dold));
    end
    
end
