function [] = SingleTest()
%% Example Test Function
% Initialize parameters for the simulation and run localization algorithm

clear variables;
clc;
close all;

%% Main simulation parameters
% Area of Interest
ROI = 100;
% Grid Element Size
gridSize = 5;
% Number of Sensors
N_s = 3;
% Transmit Power of the Emitter (Watts)
P_T = 1;
% Assumed Transmit Power (Watts)
P_E = 1;
% Shadow Spread (dB)
sigma = 0;
% Path Loss Exponent (Actual)
alpha_actual = 3.5;
% Path Loss Exponent (Assumed)
alpha_assumed = 3.5;
% Receiver Sensitivity (Watts)
recSens = -Inf;

%% Output Flags
% Display placements/estimation
dispON = 1;
% Time the Script
useTime = 0;

% TEST Purposes
% Corner Placement
% sPos(:,1) = [0.001;100];
% sPos(:,2) = [100;0.001];    
% sPos(:,3) = [100;100];
% sPos(:,4) = [0.001;0.001]; 
% Side Placement
% sPos(:,1) = [0.001;50];
% sPos(:,2) = [50;100];    
% sPos(:,3) = [50;0.001];
% sPos(:,4) = [100;50]; 
% Array Placement
% sPos(:,1) = [50;0.001];
% sPos(:,2) = [50;33];
% sPos(:,3) = [50;66];
% sPos(:,4) = [50;100];
% Triangle Placement
% sPos(:,1) = [0.001;0.001];
% sPos(:,2) = [0.001;50];
% sPos(:,3) = [100;50];
% sPos(:,4) = [0.001;100];
% ePos(:,1) = [71.01;70.01];
% Errored Placement at GridSize 9 NonLinear KE
% sPos(:,1) = [91.3337;15.2378];
% sPos(:,2) = [82.5817;53.8342];    
% sPos(:,3) = [44.2678;10.6653];
% sPos(:,4) = [96.1898;0.46342];  
% ePos(:,1) = [77.491;81.7303];
% Errored Placement at GridSize 5
% sPos(:,1) = [57.183;28.6018];
% sPos(:,2) = [69.9134;79.6258];
% sPos(:,3) = [44.1589;44.6216];
% sPos(:,4) = [46.5662;27.9039];
% ePos(:,1) = [67.5375;90.3665];
% Errored Placement at GridSize 5 for KE with distances
% sPos(:,1) = [63.6512;48.0362];
% sPos(:,3) = [61.7058;65.3223];
% sPos(:,4) = [54.9829;61.2566];
% sPos(:,5) = [58.344;82.902];
% ePos(:,1) = [72.5908;9.6768];
% Random placement
% sPos(:,1) = [82.5314;8.347];
% sPos(:,2) = [13.3171;17.3389];
% sPos(:,3) = [39.0938;83.138];
% sPos(:,4) = [39.9258;52.6876];
% ePos(:,1) = [41.6799;65.686];
% Goldoni T1 placement
sPos(:,1) = [17;2];
sPos(:,2) = [9.5;10];
sPos(:,3) = [2;2];
ePos(:,1) = [17;9.5];

% Assign Sensor Locations
%sPos = 0;
assignS = sPos;

% Assign Emitter Locations
%ePos = 0;
assignE = ePos;                  % Assign Sensors, 0 otherwise

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
% FGS_SE( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual, alpha_assumed...
%     , recSens ...
%     , plotON, dispON, tableON...
%     , useTime, useKE...
%     , assignS, assignE, centerSensors, centerEmitters);    
Trilateration( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
    , alpha_assumed, recSens ...
    , dispON ...
    , useTime, assignS, assignE);  
MinMax( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
    , alpha_assumed, recSens ...
    , dispON ...
    , useTime, assignS, assignE);
MaximumLikelihood( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
    , alpha_assumed, recSens ...
    , dispON ...
    , useTime, assignS, assignE);
end