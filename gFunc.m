function c = gFunc(tetha,pV,vRV,r,Positions)
% evaluates function g in equation (15) in 
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6475197
    NK = numel(r);
    h = sqrt(sum((Positions - repmat(tetha,1,NK)).^2))';
    v = r - h;
    [X,Y] = meshgrid(v,vRV(1:end-1));
    [~,ind] = min(abs(X-Y));
    c = -sum(log(pV(ind)));
end