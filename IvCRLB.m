% Copyright 2017 Alireza Zaeemzadeh
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Please cite the following reference if you use this code
% Robust Target Localization Based on Squared Range Iterative Reweighted Least Squares
% Alireza Zaeemzadeh, Mohsen Joneidi, Behzad Shahrasbi and Nazanin Rahnavard 
% 14th IEEE International Conference on Mobile Ad hoc and Sensor Systems (MASS) 2017
%
% Please report any bug at zaeemzadeh -at- knights -dot- ucf -dot- edu 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Iv = IvCRLB(noiseSTD,beta,outliertype,outlierparam,Nc,Gridsize)
% noiseSTD : Gaussian noise standard deviation
% beta: contamination ratio 0 < beta < 1
% outliertype: 'Gaussian' or 'Uniform'
% outlierparam: parameter of the outlier distribution
% Nc : Samples generated from the distribution

% p_v(v)  = (1-beta)*Gaussian(0,noiseSTD) + beta*Outlier(outlierparam)
    res = Gridsize/1000;            % resolution of the distribution
    v = -sqrt(2)*Gridsize:res:sqrt(2)*Gridsize;       % support of the distribution
    Gaussian = normpdf(v,0,noiseSTD);
    if strcmp(outliertype,'Uniform')
        Outlier = pdf('Uniform',v,outlierparam(1),outlierparam(2));
    elseif strcmp(outliertype,'Gaussian')
        Outlier = pdf('Normal',v,outlierparam(1),outlierparam(2));        
    end
    pv = (1-beta)*Gaussian + beta*Outlier;
    gradpv = diff(pv)/res;
    
    pvCDF = cumsum(pv)*res;
    Iv = 0;
    for n = 1:Nc
        % draw a sample from pv
        p = rand;
        [~,ind] = min(abs(p-pvCDF));
        
        Iv = Iv + (gradpv(ind)/pv(ind))^2;
    end
    Iv = Iv/Nc;
%     figure(4);plot(v,pv);
end
% reference:
% TOA-Based Robust Wireless Geolocation and
% Cramér-Rao Lower Bound Analysis in Harsh
% LOS/NLOS Environments
% http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6475197&tag=1