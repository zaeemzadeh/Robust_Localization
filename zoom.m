function alpha = zoom(alphalo,alphahi,s,tetha,pV,vRV,r,Positions)
% zoom algorithm 
% J. Nocedal and S. J. Wright  Numerical Optimization,  1999 :Springer-Verlag
% Algorithm 3.6
    c1 = 0;
    c2 = 1;
    maxIter = 20;
    for i = 1:maxIter
        alpha = (alphalo+alphahi)/2;

        Phi = gFunc(tetha+alpha*s,pV,vRV,r,Positions);

        Phi0 = gFunc(tetha,pV,vRV,r,Positions) ;

        deltaalpha = .1;
        PhiDiff0 = gFunc(tetha + deltaalpha*alpha*s,pV,vRV,r,Positions) - gFunc(tetha,pV,vRV,r,Positions);
        PhiDiff0 = PhiDiff0/(deltaalpha*alpha);

        if (Phi > Phi0 + c1*alpha*PhiDiff0) || Phi > gFunc(tetha+alphalo*s,pV,vRV,r,Positions)
            alphahi = alpha;
        else
            PhiDiff = gFunc(tetha + (1+deltaalpha)*alpha*s,pV,vRV,r,Positions) - gFunc(tetha,pV,vRV,r,Positions);
            PhiDiff = PhiDiff/(deltaalpha*alpha);
            if abs(PhiDiff) <= -c2*PhiDiff0
                break;
            end
            if PhiDiff*(-alphalo+alphahi)>=0
                alphahi = alphalo;
            end
            alphalo = alpha;
        end
    end
end