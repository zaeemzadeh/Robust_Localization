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
function CRLB = CRLBpos(theta,X,Y,Iv)
% Cramer Rao Lower Bound for
% Robust Target Localization Based on Squared Range Iterative Reweighted Least Squares
% Alireza Zaeemzadeh, Mohsen Joneidi, Behzad Shahrasbi and Nazanin Rahnavard 
% 14th IEEE International Conference on Mobile Ad hoc and Sensor Systems (MASS) 2017
% theta : target position
% X and Y : sensor positions
% Iv : intrinsic accuracy
% CRLB : Cramer-Rao Lower Bound

    htheta = sqrt(sum([X - theta(1) Y - theta(2)].^2,2));
    Iv = Iv * eye(numel(X));

    H = [(theta(1) - X)./htheta (theta(2) - Y)./htheta]';
    
    F = H * Iv * H';                % Fisher information matrix
    
    CRLB = sqrt(trace(inv(F)));
end

% reference:
% TOA-Based Robust Wireless Geolocation and
% Cramér-Rao Lower Bound Analysis in Harsh
% LOS/NLOS Environments
% http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6475197&tag=1