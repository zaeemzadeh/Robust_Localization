function [pV,vRV] = KDE(v,noiseSTD)
% Construct an estimate of the true PDF via the nonparametric AKDE
%  http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6475197 (Algorithm III)
    NK = numel(v);
    % find a pilot density estimator
    w0 = 0.79*iqr(v)*(NK^(-0.2));
    vRV = linspace(-10*noiseSTD,10*noiseSTD,100);
    p0 = zeros(1,numel(vRV));
    for m = 1:NK
        p0 = p0 + exp( -0.5*((vRV-v(m))/w0).^2 ) /(sqrt(2*pi)*w0);
    end
    p0 = p0/NK;
    p0 = p0/max(p0);
    p0(p0<1e-3) = 1e-3;
    p0 = p0/sum(p0);

    % define local bandwidth
    beta = 0.5;
    [X,Y] = meshgrid(v,vRV);
    [~,ind] = min(abs(X-Y));
    lambda = (p0(ind) / geomean(p0) ).^beta;
    % global bandwidth is determined using the least-squares cross-validation technique [33].
    w = 0.1*w0:0.1*w0:10*w0;
    M1 = zeros(1,numel(w));
    for i = 1:numel(w)
        XiXj = repmat(v,1,NK) - repmat(v,1,NK)';
        Kstar = exp( -0.5*(XiXj./2*w(i)).^2 )/(sqrt(2*pi)*2) - 2* exp( -0.5*(XiXj./w(i)).^2 )/(sqrt(2*pi));
        M1(i) = sum(Kstar(:))/(w(i)*NK^2) + 2/(sqrt(2*pi)*NK*w(i));
    end
    [~,ind] = min(M1);
    w = w(ind);
    % construct the adaptive kernel density estimator
    pV = zeros(1,numel(vRV));
    for m = 1:NK
        pV = pV + exp( -0.5*((vRV-v(m))./(w*lambda(m))).^2 ) /(sqrt(2*pi)*w*lambda(m));
    end
    pV = pV/NK;
end