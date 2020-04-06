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
%%
clear all
rng('default')
%%
sweep_parameter = 10:10:60;    %parameter to sweep ****************
RandomRuns = 50;       % ********************

ErrorSRLS = zeros(size(sweep_parameter));
ErrorSRL0 = zeros(size(sweep_parameter));
ErrorSRBCD = zeros(size(sweep_parameter));
ErrorSRHybrid = zeros(size(sweep_parameter));
ErrorLowerBound = zeros(size(sweep_parameter));

BiasSRLS = zeros(2,numel(sweep_parameter));
BiasSRL0 = zeros(2,numel(sweep_parameter));
BiasSRBCD = zeros(2,numel(sweep_parameter));
BiasSRHybrid = zeros(2,numel(sweep_parameter));

TimeSRLS = zeros(size(sweep_parameter));
TimeSRL0 = zeros(size(sweep_parameter));
TimeSRBCD = zeros(size(sweep_parameter));
TimeSRHybrid = zeros(size(sweep_parameter));

NumItersSRL0 = zeros(size(sweep_parameter));
NumItersSRBCD = zeros(size(sweep_parameter));
NumItersSRHybrid = zeros(size(sweep_parameter));
%%
for p = 1:numel(sweep_parameter)
    disp(['parameter #' num2str(p) '/' num2str(numel(sweep_parameter))]);
    %% Initialization
    error1 = zeros(1,RandomRuns);
    error3 = zeros(1,RandomRuns);
    error6 = zeros(1,RandomRuns);
    error7 = zeros(1,RandomRuns);
    CRLB = zeros(1,RandomRuns);
    
    bias1 = zeros(2,RandomRuns);
    bias3 = zeros(2,RandomRuns);
    bias6 = zeros(2,RandomRuns);
    bias7 = zeros(2,RandomRuns);
    
    time1 = zeros(1,RandomRuns);
    time3 = zeros(1,RandomRuns);
    time6 = zeros(1,RandomRuns);
    time7 = zeros(1,RandomRuns);   
    
    iters3 = zeros(1,RandomRuns);
    iters6 = zeros(1,RandomRuns);
    iters7 = zeros(1,RandomRuns);
    
    
    for seed = 1:RandomRuns
        %% Environment Parameters
        
        disp(['... Monte Carlo simulation #' num2str(seed) '/' num2str(RandomRuns)]);
        rng(seed)
        n = 2;                              % space dimension
        NumOfSensors = sweep_parameter(p); %;      % ********************
        NumofTimeSamples = 1;
        NumOfMeas = NumOfSensors*NumofTimeSamples;
        beta = 0.4;                         % contamination ratio
        NumOfOutliers = round(beta*NumOfMeas);   % ********************
        Gridsize = 4000;                      % in meters for example

        MeasNoiseSTD = 55;                   %Measurement Noise Standard Deviation in meters % ********************
        OutlierType = 'Uniform';

        SensorPosition = rand(NumOfSensors,n) * Gridsize - Gridsize/2;% ********************
%       SensorPosition = [ 100 -1300; -550 -900; 1650 -50; -1300 50; 1550 800; 350 800; -1550 800; -800 1200];

%        TargetPosition = (rand(1,n) - 0.5) * 1000 ;       % ********************
        TargetPosition = rand(1,n) * Gridsize - Gridsize/2;
%       TargetPosition = [250 250];

%         scatter(SensorPosition(:,1),SensorPosition(:,2))
    %% Outlier Parameters
        if strcmp(OutlierType,'Uniform')
            OutlierParam = [-sqrt(2)*Gridsize sqrt(2)*Gridsize]; %Outlier Measurement Range (Uniform Distribution)
        elseif strcmp(OutlierType,'Gaussian')
            OutlierParam = [380 120 ]; %Outlier Measurement Mean and STD (Gaussian Distribution)
        end
        
        if NumOfOutliers
            Outliers = randi(NumOfMeas,1,NumOfOutliers);  %Outlier Sensors
        else
            Outliers = double.empty(0,1);
        end
        
        if seed == 1
            Iv = IvCRLB(MeasNoiseSTD,beta,OutlierType,OutlierParam,1000,Gridsize);
        end
        %% Generating Measurements
        Ranges = sum((SensorPosition - repmat(TargetPosition,NumOfSensors,1)).^2,2);
        Ranges = kron(sqrt(Ranges),ones(NumofTimeSamples,1));
        SensorPosition = kron(SensorPosition,ones(NumofTimeSamples,1));
        
        % Noise model % ********************
        RangeMeas = Ranges + randn(NumOfMeas,1)*MeasNoiseSTD;  
         if strcmp(OutlierType,'Uniform')
            RangeMeas(Outliers) = Ranges(Outliers) + rand(NumOfOutliers,1)*(OutlierParam(2) - OutlierParam(1)) + OutlierParam(1);  % uniform outlier model   
        elseif strcmp(OutlierType,'Gaussian')
             RangeMeas(Outliers) = Ranges(Outliers) + randn(NumOfOutliers,1)*OutlierParam(2) + OutlierParam(1);                    % Gaussian mixture model
         end
        RangeMeas(RangeMeas<=0) = 1e-5;
        RangeMeas(RangeMeas>=sqrt(2)*Gridsize) = sqrt(2)*Gridsize;
        %% SR-LS 
        tic
        TargetSRLS = PUpositionSRLS(SensorPosition(:,1),SensorPosition(:,2),RangeMeas);
        time1(seed) = toc;
        error1(seed) = rms(TargetSRLS-TargetPosition');
        bias1(:,seed) = TargetSRLS-TargetPosition';
        

        %% SR-IRLS 
        tic
        [TargetSRL0,iters] = PUpositionSRL0(SensorPosition(:,1),SensorPosition(:,2),RangeMeas,MeasNoiseSTD);
        time3(seed) = toc;
        error3(seed) = rms(TargetSRL0-TargetPosition');
        bias3(:,seed) = TargetSRL0-TargetPosition';    
        iters3(seed) = iters;    
        %% SR-GD
        tic
        [TargetSRBCD,iters] = PUpositionSRBCD(SensorPosition(:,1),SensorPosition(:,2),RangeMeas,MeasNoiseSTD,[0;0]);
        time6(seed) = toc;
          error6(seed) = rms(TargetSRBCD-TargetPosition');
        bias6(:,seed) = TargetSRBCD-TargetPosition';        
        iters6(seed) = iters;
        %% SR-Hybrid
        tic
        [TargetHybrid,iters] = PUpositionSRHybrid(SensorPosition(:,1),SensorPosition(:,2),RangeMeas,MeasNoiseSTD);
        time7(seed) = toc;
         error7(seed) = rms(TargetHybrid-TargetPosition');
        bias7(:,seed) = TargetHybrid-TargetPosition';       
        iters7(seed) = iters; 
         %% CRLB
         CRLB(seed) = CRLBpos(TargetPosition,SensorPosition(:,1),SensorPosition(:,2),Iv);
    
    end
    ErrorSRLS(p) = mean(error1(isfinite(error1)));
    ErrorSRL0(p) = mean(error3(isfinite(error3)));
    ErrorSRBCD(p) = mean(error6(isfinite(error6)));
    ErrorSRHybrid(p) = mean(error7(isfinite(error7)));
    ErrorLowerBound(p) = mean(CRLB(isfinite(CRLB)));
    
    BiasSRLS(:,p) = mean(bias1,2);
    BiasSRL0(:,p) = mean(bias3,2);
    BiasSRBCD(:,p) = mean(bias6,2);
    BiasSRHybrid(:,p) = mean(bias7,2);
    
    TimeSRLS(p) = mean(time1);
    TimeSRL0(p) = mean(time3);
    TimeSRBCD(p) = mean(time6);
    TimeSRHybrid(p) = mean(time7);
    
    NumItersSRL0(p) = mean(iters3);
    NumItersSRBCD(p) = mean(iters6);
    NumItersSRHybrid(p) = mean(iters7);
end
%% Plot Data
figure(1);
clf;
plot(sweep_parameter,ErrorSRL0,'-og','LineWidth',4,'DisplayName',' SR-IRLS')
hold on
plot(sweep_parameter,ErrorSRLS,'-.>b','LineWidth',4,'DisplayName',' SR-LS')
plot(sweep_parameter,ErrorSRBCD,'-*c','LineWidth',4,'DisplayName',' SR-GD')
plot(sweep_parameter,ErrorSRHybrid,'--<m','LineWidth',4,'DisplayName',' SR-Hybrid')
plot(sweep_parameter,ErrorLowerBound,'--k','LineWidth',2,'DisplayName',' CRLB')

legend_handle  = legend('show');
set(legend_handle,'Interpreter','latex')
set(legend_handle,'FontSize',14)
set(legend_handle,'Location','Best')

yhandle = ylabel('RMSE');
set(yhandle,'Interpreter','latex','FontSize',16)

xhandle = xlabel('Num of Sensors');
set(xhandle,'Interpreter','latex','FontSize',16)

set(gca,'FontSize',12,'LineWidth',1)
grid on
saveas(gcf, '../results/RMSE.png');
%%
figure(2);
clf;
semilogy(sweep_parameter,TimeSRL0,'-og','LineWidth',4,'DisplayName','SR-IRLS')
hold on
semilogy(sweep_parameter,TimeSRLS,'-.^y','LineWidth',4,'DisplayName','SR-LS')
semilogy(sweep_parameter,TimeSRBCD,'-*c','LineWidth',4,'DisplayName','SR-GD')
semilogy(sweep_parameter,TimeSRHybrid,'-<m','LineWidth',4,'DisplayName','SR-Hybrid')

legend_handle  = legend('show');
set(legend_handle,'Interpreter','latex','FontSize',14,'Location','Best')

yhandle = ylabel('Time (s)');
set(yhandle,'Interpreter','latex','FontSize',16)

xhandle = xlabel('Num of Sensors');
set(xhandle,'Interpreter','latex','FontSize',16)

set(gca,'FontSize',12,'LineWidth',1)
grid on

saveas(gcf, '../results/Time.png');
