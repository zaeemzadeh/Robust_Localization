function pos = PUpositionRobustNonPara(X,Y,r,noiseSTD)
%% computes the location of the target by using Robust MDS algorithm
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6475197
% Algorithm I
% TOA-Based Robust Wireless Geolocation and Cramér-Rao Lower Bound Analysis in Harsh LOS/NLOS Environments
    %% remove faulty data
     constraint = (r < 0 )+(~isfinite(X))+ (~isfinite(Y)) ;
    X = X(constraint==0);
    Y = Y(constraint==0);
    r = r(constraint==0);
    %% Initializing algorithm parameters
    positions = [X';Y'];        
    NK = numel(r);    
    S = [-2*positions' ones(NK,1)];
    rhat = r.^2-sum(positions.^2)';
    
    tetha = (S'*S)\S'*rhat;
    tetha = tetha(1:2);
    maxIter = 20;
    ConvDelta = 1e-1;
    
    for j = 1:maxIter
        % determine residual vector
        h = sqrt(sum((positions - repmat(tetha,1,NK)).^2))';
        v = r - h;
        
        % Construct an estimate of the true PDF via the nonparametric AKDE
        [pV,vRV] = KDE(v,noiseSTD);
        
        % find the estimate by minimizing the log-likelihood function 
        tethahat = BFGS(positions,r,pV,vRV,tetha);
        if norm(tethahat-tetha) < ConvDelta
            break;
        end
        tetha = tethahat;
    end
    pos = tethahat;
end