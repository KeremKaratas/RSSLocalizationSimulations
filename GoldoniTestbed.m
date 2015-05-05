function [] = GoldoniTestbed()
%% Goldoni Testbed
% Initialize parameters for the simulation and run localization algorithms
% with the placements done in Goldoni's article. Plot a scatterplot with 25
% estimations and actual placements

clear variables;
clc;
close all;

%% Main simulation parameters
% Area of Interest
ROI = 19.5;
% Grid Element Size
gridSize = 0;
% Number of Sensors
N_s = 3;
% Transmit Power of the Emitter (Watts)
P_T = 1;
% Assumed Transmit Power (Watts)
P_E = 1;
% Shadow Spread (dB)
sigma = 2;
% Path Loss Exponent (Actual)
alpha_actual = 3.5;
% Path Loss Exponent (Assumed)
alpha_assumed = 3.5;
% Receiver Sensitivity (Watts)
recSens = -Inf;

%% Output Flags
% Display placements/estimation
dispON = 0;
% Time the Script
useTime = 0;

%% Simulation specific variables
numOfEst = 25;
% Goldoni placement
sPos(:,1) = [17;2];
sPos(:,2) = [9.5;10];
sPos(:,3) = [2;2];

% Assign Sensor Locations
%sPos = 0;
assignS = sPos;

% T1
ePosT1(:,1) = [17;9.5];
% T2
ePosT2(:,1) = [9.5;6];
% T3
ePosT3(:,1) = [9.5;9];

%% Display Simulation Parameters
% Display simulation paramaters if flag set
if(dispON)
    disp(['Region of Interest is: ' ...
            , num2str(ROI), 'x',  num2str(ROI)]);
    disp(['Region is divided into grids of size: '...
            , num2str(gridSize), 'x',  num2str(gridSize)]); 
    disp(['Number of sensors: ', num2str(N_s)]);        
    disp(['Path Loss Exponent (Actual): ', num2str(alpha_actual)]);
    disp(['Path Loss Exponent (Assumed): ', num2str(alpha_assumed)]);
    disp(['Transmit Power (Actual): ', num2str(P_T)]);
    disp(['Transmit Power (Assumed): ', num2str(P_E)]);
    disp(['Shadow Spread (dB): ', num2str(sigma)]);
end

%% Main Function Call    
estELocT1 = zeros(2,numOfEst);
estELocT2 = zeros(2,numOfEst);
estELocT3 = zeros(2,numOfEst);

estErrT1 = zeros(1,numOfEst);
estErrT2 = zeros(1,numOfEst);
estErrT3 = zeros(1,numOfEst);

for k=1:numOfEst
    [estErrT1, ~, ~, estELocT1(:,k)] ...
    = MinMax( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
                    , alpha_assumed, recSens ...
                    , dispON ...
                    , useTime, assignS, ePosT1);     
end
disp(['Mean estimated error for T1: ' num2str(mean(estErrT1))])

for k=1:numOfEst
    [estErrT2, ~, ~, estELocT2(:,k)] ...
    = MinMax( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
                    , alpha_assumed, recSens ...
                    , dispON ...
                    , useTime, assignS, ePosT2);  
end
disp(['Mean estimated error for T2: ' num2str(mean(estErrT2))])

for k=1:numOfEst
    [estErrT3, ~, ~, estELocT3(:,k)] ...
    = MinMax( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
                    , alpha_assumed, recSens ...
                    , dispON ...
                    , useTime, assignS, ePosT3);  
end
disp(['Mean estimated error for T3: ' num2str(mean(estErrT3))])

% Trilateration( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
%     , alpha_assumed, recSens ...
%     , dispON ...
%     , useTime, assignS, assignE);
% MaximumLikelihood( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
%     , alpha_assumed, recSens ...
%     , dispON ...
%     , useTime, assignS, assignE);

%% Plot known sensors and actual and estimated emitter locations
plotSELocations(ROI, gridSize, sPos, ePosT1, ePosT2, ePosT3 ...
                            , estELocT1, estELocT2, estELocT3);

end

function [] = plotSELocations(ROI, gridSize, sPos, ePosT1, ePosT2, ePosT3 ...
                            , estELocT1, estELocT2, estELocT3)
%% PLOTSELOCATIONS plots sensors and emittor on figure to show locations
% Required inputs are;
% - Region of interest axis length (ROI)
% - Grid axis length (gridSize)
% - Sensor positions given as [x;y] (sPos)
% - Emitter positions given as [x;y] (ePos)

figure();

% Create whole area
hold on;
rectangle('Position', [0, 0, 19.5, 12]);
% Place actual sensors
scatter(sPos(1,:),sPos(2,:), 70, 'c'...
            , 'MarkerFaceColor', [1 0 0] ...
            , 'MarkerEdgeColor', [1 0 0]);

% Place actual emitters also determine color for each
c(1,:)= [round(rand), round(rand), round(rand)];
while isequal(c(1,:),[1 1 1])
    c(1,:)= [round(rand), round(rand), round(rand)];
end

scatter(ePosT1(1,1),ePosT1(2,1), 70, 's' ...
            , 'MarkerFaceColor', [0 0 0]...
            , 'MarkerEdgeColor', [0 0 0]);
scatter(ePosT2(1,1),ePosT2(2,1), 70, 's' ...
            , 'MarkerFaceColor', [0 0 0]...
            , 'MarkerEdgeColor', [0 0 0]);
scatter(ePosT3(1,1),ePosT3(2,1), 70, 's' ...
            , 'MarkerFaceColor', [0 0 0]...
            , 'MarkerEdgeColor', [0 0 0]);        

% Place estimated emitter locations
for k=1:length(estELocT1(1,:))
    hT1 = scatter(estELocT1(1,k), estELocT1(2,k), 60, '*'...
                , 'MarkerFaceColor', 'none'...
                , 'MarkerEdgeColor', [1 0 0]);
end
for k=1:length(estELocT2(1,:))
    hT2 = scatter(estELocT2(1,k), estELocT2(2,k), 60, '+'...
                , 'MarkerFaceColor', 'none'...
                , 'MarkerEdgeColor', [0 0 1]);
end
for k=1:length(estELocT3(1,:))
    hT3 = scatter(estELocT3(1,k), estELocT3(2,k), 60, 'p'...
                , 'MarkerFaceColor', 'none'...
                , 'MarkerEdgeColor', [0 1 0]);
end
hold off;

% Set title, labels and legend
title('Sensor and Emitter placements for FGS_{SE}');
xlabel('x');
ylabel('y');
legend([hT1 hT2 hT3], {'T1', 'T2', 'T3'});

end