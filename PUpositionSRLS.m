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
function pos = PUpositionSRLS(X,Y,R)
% finding the position of PUs using SRLS method in
% Exact and Approximate Solutions of Source Localization Problems
% X : Sensor Positions (x-coordinate)
% Y : Sensor Positions (y-coordinate)
% R : Range Measurements
% pos : Target Position
%%
    constraint = (R < 0 )+(~isfinite(X))+ (~isfinite(Y));
    X = X(constraint==0);
    Y = Y(constraint==0);
    R = R(constraint==0);
    %% scaling the measurements to avoid badly scaled matices
    normalize = mean(abs([X(:);Y(:)]));
    X = X./normalize;
    Y = Y./normalize;
    R = R./normalize;
    %% initialization
    positions = [X';Y'];        
    m = numel(R);
    A = [-2*positions' ones(m,1)];
    b = R.^2-sum(positions.^2)';
    D = [eye(2) zeros(2,1) ; zeros(1,2) 0];
    f = [zeros(2,1);-0.5];
    %% solving the GTRS problem
    fun = @(lambda) ((A'*A + lambda*D)\(A'*b-lambda*f))'*D*((A'*A + lambda*D)\(A'*b-lambda*f)) + 2*f'*((A'*A + lambda*D)\(A'*b-lambda*f)); % function
%     options = optimset('Display','iter'); % show iterations
    lambdahat = fzero(fun,0);
    lambdal = max(-diag(A'*A)./diag(D));          % http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.31.2504&rep=rep1&type=pdf PAGE14
    if (lambdahat<lambdal)                          % http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.31.2504&rep=rep1&type=pdf
        lambdahat = lambdal;
    end
    yhat = (A'*A + lambdahat*D)\(A'*b-lambdahat*f);

    pos = yhat(1:2);
    pos = pos*normalize;
end

% references:
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.31.2504&rep=rep1&type=pdf
% http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4472183&tag=1  
% http://arxiv.org/pdf/1410.1386v2.pdf