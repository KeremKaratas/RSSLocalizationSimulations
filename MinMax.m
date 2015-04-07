function [estErr, sPos, ePos] = ...
MinMax( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual ...
    , alpha_assumed, recSens ...
    , dispON ...
    , useTime, assignS, assignE)
%% MINMAX
% This function implements Minmax method to localize the emitter. If
% sensor and emitter locations are not provided in assignS and assignE they
% are randomly placed over the Region of Interest (ROI). This random
% placement makes sure that sensors and emitters are placed in seperate
% grids.
%
% ROI - Region of interest
% gridSize - Length of the grids one side.
% N_s - Number of Sensors
% P_T - Actual transmit power of the emitter
% P_E - Assumed transmit power of the emitter
% sigma - shadowing effect in dB
% alpha_actual - Actual value of the Path Loss Exponent
% alpha_assumed - Assumed value of the Path Loss Exponent
% recSens - Receiver Sensitivity
% dispON - Flag for displaying simulation information (Set 1 to display)
% useTime - Flag for timing parts of simulation (Set 1 to time)
% assignS - Predefined locations for sensors (0 to randomly assign)
% assignE - Predefined locations for emitter (0 to randomly assign)
                
if useTime
    ccsTime = tic();
end

%% Placing Sensors/Emitter and Creating Sensor Distance Matrix
% Making use of the vectorized coordinate system we will put the sensors
% and the emitter in the simulation and create the distance matrices.
% Sensors and the emitter is placed randomly. They are selected so that
% none of them are on the same coordinates.
%
% The number of sensors is defined here. To change the number of sensors we
% can increase N_s. Sensor and emitter locations are assigned randomly.
% They are never in the same coordinates.

if assignS
    sPos = assignS;
    
    if (dispON>0)
        for i=1:N_s        
            disp(['Sensor '...              % Display selected coordinates
                    , num2str(i), ' is located at (', num2str(sPos(1,i))...
                    ,',', num2str(sPos(2,i)), ').']);
        end
    end
else    
    % Sensor positions on (x,y) plane
    sPos = zeros(2,N_s);
    for i=1:N_s
        isDistinct = 0;                     % Boolean flag
        while ~isDistinct        
            sPos(:,i) = [rand()*ROI;rand()*ROI];                     
            isDistinct = 1;                 % True unless proven otherwise
            for j=1:i-1
                if (abs(sPos(1,i)-sPos(1,j))<gridSize ...
                    && abs(sPos(2,i)-sPos(2,j))<gridSize)
                    isDistinct = 0;         % Shown that it is not distinct
                end                
            end
        end
        if (dispON>0)
            disp(['Sensor '...              % Display selected coordinates
                    , num2str(i), ' is located at (', num2str(sPos(1,i))...
                    ,',', num2str(sPos(2,i)), ').']);
        end
    end
end
if assignE
    ePos = assignE;
    if (dispON>0)
        disp(['Emitter is placed at ('...   % Display selected coordinates
                 , num2str(ePos(1,1)) ,',', num2str(ePos(2,1)), ').']);
    end
else    
    % Emitter position on (x,y) plane
    ePos = zeros(2,1);
    isDistinct = 0;                         % Boolean flag
    while ~isDistinct        
        ePos(:,1) = [rand()*ROI;rand()*ROI];       
        isDistinct = 1;                     % Unless proven otherwise, true
        for i=1:N_s
            if (abs(sPos(1,i)-ePos(1,1))<gridSize ...
                && abs(sPos(2,i)-ePos(2,1))<gridSize)
                % Shown that it is not distinct
                isDistinct = 0;
            end
        end
    end

    if (dispON>0)
        disp(['Emitter is placed at ('...   % Display selected coordinates
                 , num2str(ePos(1,1)) ,',', num2str(ePos(2,1)), ').']);
    end
end

if useTime
    disp(['Placing sensors and emitter took '...
                , num2str(toc(ccsTime)), ' seconds']);
    crssTime = tic();
end

%% Creating Received Signal Strength Values
% These variables are required to create the signal model 
% $m_i =  (P_T) (d_i)^{- \alpha}$ in Watts where 
% $1 \geq i \geq N_s$, $$P_T$ is the the transmit power of the emitter, 
% $\alpha$ is the path loss exponent and 
% $d_i=\sqrt{(x_i-x_0)^2 + (y_i-y_0)^2}$ is the distance between the
% transmitter and the ith sensor.
%
% The sensor's experience log-normal shadowing. If the fast fading effects
% are sufficiently averaged over time then the resulting unknown measured
% power from the emitter to the ith sensor is given as 
% $r_i = 10^{log_{10} (m_i) + \frac{\omega_i}{10}}$ where
% $\omega_i$ is a normal random variable with
% mean of 0 and a variance of $\sigma^2$.

% Sensor and emitter Distance 
d = sqrt((sPos(1,:)-ePos(1,:)).^2 + (sPos(2,:)-ePos(2,:)).^2);

% Received signal in Watts
m = P_T * d.^(-alpha_actual);
r = m;
for i=1:N_s
    r(i) = 10^(log10(m(i)) + normrnd(0,sigma)/10);
    % If received signal strength value is smaller than receiver
    % sensitivity we assume we didn't receive it    
    if(r(i)<recSens)
        r(i) = recSens;   % TODO: Change here to properly implement recSens
        if(dispON)
            disp([num2str(i), 'th Sensor is ignored.'])
        end
    end
end

if useTime
    disp(['Creating RSS took ' num2str(toc(crssTime)) ' seconds']);
    ctTime = tic();
end

%% Emitter Location Estimation
% Distances are estimated from each sensor. A square is formed using sensor
% location as the center and sides twice the length of the estimated
% distance. Estimation of the emitter is done as the center of the
% quadrilateral formed by the intersections of the squares.

% Estimate the distance from sensors to the emitter using FSPL
estDist = (r./P_E).^(-1/alpha_assumed);

% Create a square by finding the intersection of circles that centers known
% sensor locations with its related estimated distance.
xLine = [max(sPos(1,:)-estDist) min(sPos(1,:)+estDist)];
yLine = [max(sPos(2,:)-estDist) min(sPos(2,:)+estDist)];

% Find the center of the square. Estimated location of emitter is here.
ex = mean(xLine);
ey = mean(yLine);

if dispON>0
    if length(ex)>1
        disp([num2str(length(ex)), ' possible solutions are found']);
    end
end
if dispON>0
    for i=1:length(ex)
        disp(['Emitter found at (', num2str(ex(i))...
                , ',', num2str(ey(i)), ').']);
    end
end
ex = ex(1);
ey = ey(1);

% Assign estimated emitter location
estELoc = [ex;ey];

% Now we have both the emitter location and the estimation, we can find
% the error.
estErr = sqrt((ePos(1,:)-ex).^2 + (ePos(2,:)-ey).^2);
if dispON>0
    disp(['Estimation error is ', num2str(estErr)  , ' meters.']);
end

if(dispON)
    plotSELocations(ROI, gridSize, sPos, ePos, estELoc, estDist);
end
                            
if useTime
    disp(['Calculating theta took ' num2str(toc(ctTime)) ' seconds']);
    disp(['Total elapsed time is ' num2str(toc(ccsTime)) ' seconds']);
end

end

function [] = plotSELocations(ROI, gridSize, sPos, ePos, estELoc, estDist)
%% PLOTSELOCATIONS plots sensors and emittors on figure to show locations
% Required inputs are;
% - Region of interest axis length (ROI)
% - Grid axis length (gridSize)
% - Sensor positions given as [x;y] (sPos)
% - Emitter positions given as [x;y] (ePos)

figure();

% Create whole area
hold on;

% Place actual sensors
scatter(sPos(1,:),sPos(2,:), 70, 'c'...
            , 'MarkerFaceColor', [1 0 0] ...
            , 'MarkerEdgeColor', [1 0 0]);

% Place actual emitters also determine color for each
c(1,:)= [round(rand), round(rand), round(rand)];
while isequal(c(1,:),[1 1 1])
    c(1,:)= [round(rand), round(rand), round(rand)];
end

scatter(ePos(1,1),ePos(2,1), 70, 's' ...
            , 'MarkerFaceColor', [0 0 0]...
            , 'MarkerEdgeColor', [0 0 0]);

% Place estimated emitter locations
scatter(estELoc(1,1), estELoc(2,1), 70, 'd'...
            , 'MarkerFaceColor', [0 0 0]...
            , 'MarkerEdgeColor', [0 0 0]);

% plot the squares and the point found
t = linspace(0,2*pi)';
for i = 1:length(sPos(1,:))
    rectangle('Position', [sPos(1,i)-estDist(i) sPos(2,i)-estDist(i) ...
                            estDist(i)*2 estDist(i)*2]);
end        
        
hold off;
% Set grid on with the given size
% tickValues = 0:gridSize:ROI+1;
% set(gca,'YTick',tickValues)
% set(gca,'XTick',tickValues)
% grid on;

% Set title, labels and legend
title('Sensor and Emitter placements for FGS_{SE}');
xlabel('x');
ylabel('y');
legend('Sensors', 'Emitter', 'Estimated Emitter');

end