function g = gFuncGradient(tetha,pV,vRV,r,Positions)
% evaluates gradient of the function g in equation (15) in 
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6475197
    NK = numel(r);
    h = sqrt(sum((Positions - repmat(tetha,1,NK)).^2))';
    v = r - h;
    pVdiff = diff(pV)./diff(vRV);

    [X,Y] = meshgrid(v,vRV(1:end-1));
    [~,ind] = min(abs(X-Y));

    g = (pVdiff(ind).^2)./pV(ind);
    g = repmat(g,2,1) .* (repmat(tetha,1,NK) - Positions) ./ repmat(h',2,1) ;
    g = sum(g,2);       %gradient 
end