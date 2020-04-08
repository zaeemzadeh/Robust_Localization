function tetha = BFGS(Positions,r,pV,vRV,tetha)
% BFGS quasi-Newton method with a cubic line search
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6475197
% Algorithm 2
    H = eye(2);
    maxIter = 20;
    for l = 1:maxIter
        % obtain search diraction
        g = gFuncGradient(tetha,pV,vRV,r,Positions);       %gradient   
        s = - H*g;      % descent direction
        s = s./norm(s);
        
        % find step size  (linesearch algorithm)
        alpha0 = 0;      %stepsize
        alphamax = 100;
        for i = 1:maxIter
            alpha = rand*(alphamax-alpha0) + alpha0;
            %evaluate Phi(alpha) = g(tetha_i +alpha*s_i )
            c1 = 1e-4;
            c2 = 1e-3;
            Phi = gFunc(tetha + alpha*s,pV,vRV,r,Positions);
            % evaluate Phi'(0)
            deltaalpha = .1;
            PhiDiff0 = gFunc(tetha + deltaalpha*alpha*s,pV,vRV,r,Positions) - gFunc(tetha,pV,vRV,r,Positions);
            PhiDiff0 = PhiDiff0/(deltaalpha*alpha);
            
            Phi0 = gFunc(tetha,pV,vRV,r,Positions) ;
            
            if ((Phi > Phi0 + c1*alpha*PhiDiff0) || (Phi >= gFunc(tetha+alpha0*s,pV,vRV,r,Positions) && i>1))
                alpha = zoom(alpha0,alpha,s,tetha,pV,vRV,r,Positions);
                break;
            end
            % evaluate Phi'(alpha)
            deltaalpha = 0.1;
            PhiDiff = gFunc(tetha + (1 + deltaalpha)*alpha*s,pV,vRV,r,Positions) - gFunc(tetha + alpha*s,pV,vRV,r,Positions);
            PhiDiff = PhiDiff/(deltaalpha*alpha);
            if abs(PhiDiff) <=-c2*PhiDiff0
                break;
            end
            
            if PhiDiff >= 0
                alpha = zoom(alpha,alpha0,s,tetha,pV,vRV,r,Positions);break;
            end 
        end
        
        % Update approximate Hessian Matrix
        delta = alpha*s;
        gamma = gFuncGradient(tetha+alpha*s,pV,vRV,r,Positions) - gFuncGradient(tetha,pV,vRV,r,Positions);
        H = H + (1 + (gamma'*H*gamma)/(delta'*gamma))*(delta*delta')/(delta'*gamma) - (delta*gamma'*H+H*gamma*delta')/(delta'*gamma);
        if norm(alpha*s) < 1e-6
            tetha = tetha + alpha*s;break;
        end
        % update estimate
        tetha = tetha + alpha*s;
    end
    
end