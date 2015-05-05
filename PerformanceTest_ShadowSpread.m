function [] = PerformanceTest_ShadowSpread()
%% Performance test for shadow spread
% Tests the performance of a localization technique for different shadow
% spread values.

clear variables;
clc;
close all;

%% Main simulation parameters
% Area of Interest
ROI = 1000;
% Grid Element Size
gridSize = 1;
% Number of Sensors
N_s = 50;
% Transmit Power of the Emitter (Watts)
P_T = 1;
% Assumed Transmit Power (Watts)
P_E = 1;
% Path Loss Exponent (Actual)
alpha_actual = 3.5;
% Path Loss Exponent (Assumed)
alpha_assumed = 3.5;
% Receiver Sensitivity (Watts)
recSens = -Inf;

%% Output Flags
% Number of zoomed plots
plotON = 0;
% Display placements/estimation
dispON = 0;
% Assign Sensor Locations
assignS = 0;
% Assign Emitter Locations
assignE = 0;
% Time the Script
useTime = 0;

%% Simulation Spesific variables and Console Output of parameters

% Number of Function Calls
trialSize = 1000;

% Shadow Spread Values
SS = 0:2:10;

% Calculation of mean estimated error values for different shadow spread.
disp(' ');
disp('Trilateration Mean Estimation Errors for Different Shadow Spread Values');
disp(' ');
disp('Simulation parameters are given below');
disp(['ROI = ', num2str(ROI), ', N_s = ', num2str(N_s) ...
    , ', trialSize = ', num2str(trialSize)]);
disp(['P_T = ', num2str(P_T), ', P_E = ', num2str(P_E) ...
    , ', gridSize = ', num2str(gridSize), ',']);
disp(['Path Loss Exponent (Actual - Assumed) = ' ...
    , num2str(alpha_actual), ' - ', num2str(alpha_assumed)]);
disp(' ');
disp('Mean error per shadow spread value;');

% Initialize simulation variables
estErr = zeros(trialSize,length(SS));
meanEstErrTrilateration = zeros(1,length(SS));
meanEstErrMinMax = zeros(1,length(SS));
meanEstErrMaximumLikelihood= zeros(1,length(SS));

%% Trilateration Simulation
disp(' ');
disp('Trilateration Simulation');
disp(' ');

for sidx=1:numel(SS)
    sigma=SS(sidx);
    parfor i=1:trialSize
        estErr(i,sidx)...
        = Trilateration( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
                        , alpha_assumed, recSens ...
                        , dispON ...
                        , useTime, assignS, assignE);
    end
    meanEstErrTrilateration(sidx) ...
        = sum(estErr(:,sidx))/trialSize;
    disp(['Shadow Spread ', num2str(sigma) ...
            , ' is complete with mean error: ' ...
            , num2str(meanEstErrTrilateration(sidx)), ' meters.']);
end


%% MinMax Simulation
disp(' ');
disp('MinMax Simulation');
disp(' ');

for sigma=SS    
    for i=1:trialSize            
        estErr(i,find(SS==sigma))...                
        = MinMax( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
                        , alpha_assumed, recSens ...
                        , dispON ...
                        , useTime, assignS, assignE);
    end
    meanEstErrMinMax(find(SS==sigma)) ...
        = sum(estErr(:,find(SS==sigma)))/trialSize;
    disp(['Shadow Spread ', num2str(sigma) ...
            , ' is complete with mean error: ' ...
            , num2str(meanEstErrMinMax(find(SS==sigma))), ' meters.']);            
end

%% Maximum Likelihood Simulation
disp(' ');
disp('Maximum Likelihood Simulation');
disp(' ');

for sigma=SS    
    for i=1:trialSize            
        estErr(i,find(SS==sigma))...                
        = MaximumLikelihood( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
                        , alpha_assumed, recSens ...
                        , dispON ...
                        , useTime, assignS, assignE);
    end
    meanEstErrMaximumLikelihood(find(SS==sigma)) ...
        = sum(estErr(:,find(SS==sigma)))/trialSize;
    disp(['Shadow Spread ', num2str(sigma) ...
            , ' is complete with mean error: ' ...
            , num2str(meanEstErrMaximumLikelihood(find(SS==sigma))) ...
            , ' meters.']);
end
%% Plotting the Performance
hold on;
plot(SS, meanEstErrTrilateration, 'x-');
plot(SS, meanEstErrMinMax, 'v-');
plot(SS, meanEstErrMaximumLikelihood, 'p-');
grid on;
title('Mean Estimation Error for Different Shadow Spread Values');
xlabel('Shadow spread (dB)');
ylabel('Mean Estimation error (m)');
legend('Trilateration' ...
       , 'MinMax' ...
       , 'Maximum Likelihood');
hold off;
end