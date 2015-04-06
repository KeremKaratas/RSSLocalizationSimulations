%% Iterative Grid Search for RSS-Based Emitter Localization
% Main function implements the Full Grid search method provided in the 
% article given in the title. The variables required as input to the
% function is described below. Detailed descriptions of these variables are
% given in under the "Main Function" section.
%
% <html>
% <table border=1>
% <tr><td>Name</td>       <td>Description</td></tr>
% <tr><td>ROI</td>        <td>Region of Interest</td></tr>
% <tr><td>gridSize</td>   <td>Grid Size</td></tr>
% <tr><td>N_s</td>        <td>Number of Sensors</td></tr>
% <tr><td>P_T</td>        <td>Real Emitter Transmit Power</td></tr>
% <tr><td>sigma</td>      <td>Shadow Spread</td></tr>
% <tr><td>alpha</td>      <td>Path Loss Exponent</td></tr>
% <tr><td>plotON</td>     <td>Plot Zoom Level</td></tr>
% <tr><td>dispON</td>     <td>Display placement and estimation</td></tr>
% <tr><td>assignS</td>    <td>Assign sensors, 0 otherwise</td></tr>
% <tr><td>assignE</td>    <td>Assign emitter, 0 otherwise</td></tr>
% </table>
% </html>
%
% Dependencies for Plots:
% makedatatip.m
% subtightplot.m

function [estErr, sPos, ePos, estELoc] = ...
FGS_SE( ROI, gridSize, N_s, P_T, P_E, sigma, alpha_actual, alpha_assumed...
    , recSens ...
    , plotON, dispON, tableON...
    , useTime, useKE...
    , assignS, assignE, centerSensors, centerEmitters)                  
%% Main Function
% This function implements a single estimation of an emitter using sensors
% that are placed randomly. First a gridded coordinate system is created.
% Then sensor and emitter locations are assigned. After placement received
% signal power values are created. Using these values and the known sensor
% locations emitters location is estimated. The error of estimation is
% found and returned at the end of the function.
%
% If plotON is given as 0 nothing will be displayed and there will be no
% plots. Changing it to something larger than 1 will create a plot multiple
% subplots that shows the cost function with different vertical zoom levels
% . Also sensor and emitter placements and the estimation is displayed.

                
if useTime
    ccsTime = tic();
end

%% Creating the coordinate system
% Here we create the 2D coordinate system as two matrices each including
% their associated indexes. This plane will be our Region of Interest(ROI).
% ROI is divided into grids to minimize computation.
%
% Vectorization in MATLAB is important as reallocation of memory takes a
% lot of time in simulations therefore everything should be predefined for
% best performance.

% Locations on x axis
u = 0:1:ROI-1;
% Locations on y axis
v = 0:1:ROI-1;
% 2D plane indices for x and y
[x,y] = meshgrid(u, v);

newROI = ceil(ROI/gridSize);                % ROI size with grid
xGrid = zeros(newROI ,newROI );
yGrid = zeros(newROI ,newROI );

for k=0:gridSize:ROI-1
    xGrid(:,(k/gridSize+1)) ...
        = x(1:newROI , k+1)+gridSize/2;
    yGrid((k/gridSize+1),:) ...
        = y(k+1,1:newROI)+gridSize/2;
end

% Calculate time passed from start
if useTime
    disp(['Creating the coordinate system took '...
                ,num2str(toc(ccsTime)), ' seconds']);
    pseTime = tic();
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
            if(~centerSensors)
                sPos(:,i) = [rand()*ROI;rand()*ROI];
            else
                sPos(:,i) = [xGrid(1,randi([1 newROI], 1, 1)); ...
                                yGrid(randi([1 newROI], 1, 1),1)];
            end
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
        if(~centerEmitters)
            ePos(:,1) = [rand()*ROI;rand()*ROI];
        else
            ePos(:,1) = [xGrid(1,randi([1 newROI], 1, 1)); ...
                                yGrid(randi([1 newROI], 1, 1),1)];
        end              
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

% Distance from Sensor matrices
sDist = zeros(newROI, newROI, N_s);
for i=1:N_s
    sDist(:,:,i) = sqrt((sPos(1,i)-xGrid).^2 + (sPos(2,i)-yGrid).^2);
end

if useTime
    disp(['Placing sensors and emitter took '...
                , num2str(toc(pseTime)), ' seconds']);
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

m = P_T * d.^(-alpha_actual);                  % Received signal in Watts
r = m;
for i=1:N_s
    r(i) = 10^(log10(m(i)) + normrnd(0,sigma)/10);
    % If received signal strength value is smaller than receiver
    % sensitivity we assume we didn't receive it    
    if(r(i)<recSens)
        r(i) = recSens;   
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
% Probability of observing all sensors outputs given $\theta$ is written as
% $p(r|\theta)=\prod_{i=1}^{N_s} \frac{1}{\sqrt{2 \pi \sigma^2}}
% e^{-\frac{(r_i-m_i)^2}{2 \sigma^2}}$ where $\theta = (x_0,y_0)$ is the
% emitter location parameter to be estimated. ML estimate of the emitter
% location can be obtained by maximizing the (log)likelihood function which
% requires minimization of the sum of squared differences between each
% sensor's emitter power estimate and average of all sensor estimates:
% $\theta = min_{x,y} \sum_{i = 1}^{N_s} ( 10 log_{10}(r_i (d_i)^{\alpha}))
% -\frac{1}{N_s} \sum_{i=1}^{N_s} (10 log_{10} (r_i (d_i)^{\alpha})))^2$.
% For known emitter power values the used cost function is 
% $\sum_{i = 1}^{N_s} (10 log_{10}(r_i)-10log_{10}((P_{T_E}d^{-\alpha}))^2$


innerSum =...                               % Mean of estimated power
        zeros(newROI,newROI); 
theta = zeros(newROI,newROI);       
estDist = zeros(N_s, 1);

if (useKE)        
    % $\sum_{i=1}^{N_s} ((\frac{r_i}{P_T})^\frac{1}{-\alpha} - d )^2$
    for k=1:N_s
        estDist(k) = (r(k)/P_T).^(1/-alpha_assumed);
        theta = theta + ( estDist(k) - sDist(:,:,k)).^2;
    end
    % $\sum_{i = 1}^{N_s} (10 log_{10}(r_i)-10log_{10}((P_{T_E}d^{-\alpha}))^2$
%     for i=1:N_s
%         theta = theta + (10*log10(r(i)) ...
%                 - 10*log10(P_E*sDist(:,:,i).^(-alpha_assumed))).^2;
%         
%     end
else    
    for i=1:N_s        
        innerSum = innerSum + 10*log10(r(i)*sDist(:,:,i).^alpha_assumed);
    end
    for i=1:N_s       
        theta = theta ...                      
        + (10*log10(r(i)*(sDist(:,:,i).^alpha_assumed))-(1/N_s)*innerSum).^2;            
        %+ (10*log10(r(i)*sDist(:,:,i).^alpha_assumed)).^2;            
        %+ (10*log10(r(i)*sDist(:,:,i).^alpha_assumed)-(1/(N_s-1))*innerSum).^2;
    end
end

% Finding of the indices of the minimum value of the $\theta$ function 
% is given below.

[Iy,Ix]=find(theta==min(min(theta)));

ex = xGrid(1,Ix);
ey = yGrid(Iy,1);
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

% Draw the cost function table
% Not implemented correctly yet. Assigns the cost function to the workspace
% instead.

if tableON
    ynames = xGrid(1,:);
    xnames = yGrid(:,1);    
    assignin('base', 'ynames', ynames);
    assignin('base', 'xnames', xnames);
    assignin('base', 'theta', theta);
    assignin('base', 'sDist', sDist);
    disp('To see the table open ''theta'' variable in the workspace.')
    disp(['Emitter is minimum at the indices ['...
        , num2str(Iy(1)), ','...
        , num2str(Ix(1)), '] of the theta variable.']);
end

% Plot the resulting cost function
if(plotON~=0)    
    % We will use this figure to plot the cost function 3D
    figure('name','Emitter Location Estimation','numbertitle','off');

    % Find how many subplots are neccessary
    plotSize = 0;
    while(abs(plotON)>plotSize^2)
        plotSize=plotSize+1;
    end

    h=zeros(1,plotSize);
    for i=1:abs(plotON)        
        if (exist('subtightplot', 'file'))
            subtightplot(plotSize,plotSize,i...
                            ,[0.1 0.1], [0.1 0.1], [0.1 0.1]);
        else
            subplot(plotSize,plotSize,i)
        end
        % Change large values to NaN for better visual
        visual = theta;           
        if(plotON>0)
            for j=1:i-1
                visual(visual>=(max(visual(:))-min(visual(:)))/2) = NaN;
            end
        end
        if(plotON<0)
            visual(visual>=nanmean(visual(:))) = NaN;
            for j=1:i
                visual(visual>=nanmean(visual(:))) = NaN;
            end
        end       
        colormap([jet(128);gray(128)])
        h(i) = surf(xGrid,yGrid,visual);

        % Initially, both CDatas are equal to cost.
        color_S = 128; % 128-elements is each colormap

        % CData for surface
        cmin = min(visual(:));
        cmax = max(visual(:));
        C1 = min(color_S,round((color_S-1)*(visual-cmin)/(cmax-cmin))+1); 

        % CData for pcolor
        C2 = 128+C1;

        % Update the CDatas for each object.
        set(h(i),'CData',C1);

        % Change the CLim property of axes so that it spans the 
        % CDatas of both objects.
        caxis([min(C1(:)) max(C2(:))])

        set(h(i),'EdgeColor','none');       % This setting will make sure 
                                            % that the plot is not black 
                                            % because of default wire mesh

        % Set labels and title
        title(['Vertical Zoom Level ', num2str(i)]);
        xlabel('x');
        ylabel('y');
        zlabel('Cost/Objective Function');


        % Mark minimum on the map       
        try
            makedatatip(h(i),[Iy Ix]);
        catch            
            disp('Minimum value not shown as makedatatip.m is not found.');
        end
    end

    % Enlarge figure
    set(gcf, 'units','normalized','outerposition',[0.05 0.13 0.9 0.8]);  
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
axis([0 ROI 0 ROI]);
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

% plot the circles and the point found
if(estDist)
    t = linspace(0,2*pi)';
    for i = 1:length(sPos(1,:))
        plot(sPos(1,i) + estDist(i)*cos(t),sPos(2,i) + estDist(i)*sin(t),'-')
    end            
end

hold off;
% Set grid on with the given size
tickValues = 0:gridSize:ROI+1;
set(gca,'YTick',tickValues)
set(gca,'XTick',tickValues)
grid on;

% Set title, labels and legend
title('Sensor and Emitter placements for FGS_{SE}');
xlabel('x');
ylabel('y');
legend('Sensors', 'Emitter', 'Estimated Emitter');

end

function varargout = makedatatip(hObj,index)
%MAKEDATATIP  Adds data tips to specified data points of graphical objects.
%
%  MAKEDATATIP(HOBJ,INDEX) adds a datatip to the graphical object HOBJ at
%  the data point defined by INDEX.
%
%  HOUT = MAKEDATATIP(...) returns handles to the created datatips.
%
%  If HOBJ is 1-dimensional, INDEX can be of any size and is assumed to be
%  a linear index into the data contained by HOBJ.  HOUT will be of the
%  same size as INDEX.
%
%  If HOBJ is 2-dimensional and INDEX is N-by-2, INDEX is assumed to be a
%  set of N subscripts, and HOUT will be N-by-1.  If INDEX is of any other
%  size, it is assumed to a linear index and HOUT will be the same size as
%  INDEX.  Note that if you wish to specify 2 linear indices, ensure INDEX
%  is a column vector, else it will be assumed to be a single set of
%  subscripts.
%
% Example:
%     x = 1:10;
%     y = rand(1,10);
%     hPlot = plot(x,y);
%     makedatatip(hPlot,[3 8])
%
% Example:
%     [X,Y,Z] = peaks(30);
%     hObj = surf(X,Y,Z);
%     makedatatip(hObj,[5 8; 20 12; 22 28])
%
% Example: Add a single data tip to data point (5,25)
%     [X,Y,Z] = peaks(30);
%     hObj = surf(X,Y,Z);
%     makedatatip(hObj,[5 25])
%
% Example: Add two data tips to data points (5) and (25)
%     [X,Y,Z] = peaks(30);
%     hObj = surf(X,Y,Z);
%     makedatatip(hObj,[5; 25])
%
% Example: Add two data tips to an image
%     load mandrill
%     figure
%     hObj = image(X);
%     colormap(map) 
%     makedatatip(hObj, [103 348; 270 348])

% Author: Tim Farajian
% Release: 2.0
% Release date: 6/27/2012

% Check # of inputs
narginchk(2, 2)
nargoutchk(0, 1)

if length(hObj)~=1
  error('MAKEDATATIP:InvalidSize',...
    'HOBJ must be scalar.');
end

% Ensure hObj is valid target
if ~ishandle(hObj)
  error('MAKEDATATIP:InvalidHandle',...
    'HOBJ is an invalid handle object.');
end

isImage = strcmp(get(hObj, 'Type'), 'image'); %Determine if target is image

% Read data from hObj
try
  X = get(hObj,'XData');
  Y = get(hObj,'YData');
catch ME
  % Object must have an XData and YData property to be valid
  error('MAKEDATATIP:InvalidObjectType',...
    'Objects of class ''%s'' are not a valid targets for datatips.',...
    class(handle(hObj)))
end
try
  Z = get(hObj,'ZData');
catch ME
  % Many objects do not have a ZData property.  Some will work, some will
  % not.
  isImage = true;
end
% Ensure subscripts or indices are valid values and sizes
if isempty(index)
  return
elseif ~isnumeric(index)
  error('MAKEDATATIP:InvalidDataType',...
    'Subscript indices must be of numeric data type.')
elseif any(index(:)<1) ||...
    any(fix(index(:))~=index(:)) ||...
    any(isinf(index(:)))
  error('MAKEDATATIP:InvalidIndex',...
    'Subscript indices must be positive integers.')
elseif ~isvector(index) && ~any(size(index)==2)
  error('MAKEDATATIP:InvalidIndexMatrixSize',...
    'Subscript indices must be a vector or N-by-2 matrix.')
elseif (~isImage && isvector(X)) || size(index,2)~=2
  hDatatip = zeros(size(index));
  index = index(:);
  isLinear = true;
else
  hDatatip = zeros(size(index,1),1);
  isLinear = false;
end

% Get handle to datacursor mode object
hDataCursorMgr = datacursormode(ancestor(hObj,'figure'));

% Loop through each specified data point
for n = 1:size(index,1)
  
  % Create position vector
  if isImage && isLinear
    [i j] = ind2sub([X(2) Y(2)], index(n));
    pos = [i j 1];
  elseif isImage
    pos = [index(n, 1) index(n, 2) 1];
  elseif isempty(Z)
    pos = [X(index(n)) Y(index(n))];
  elseif isLinear
    pos = [X(index(n)) Y(index(n)) Z(index(n))];
  else
    pos = [...
      X(index(n,1),index(n,2))...
      Y(index(n,1),index(n,2))...
      Z(index(n,1),index(n,2))];
  end
  
  % Create datatip
  hDatatip(n) = createDatatip(hDataCursorMgr, hObj);
  % Specify data cursor properties
  if isImage
    set(get(hDatatip(n),'DataCursor'),'DataIndex',pos,...
      'TargetPoint',pos(1:2))
  else
    set(get(hDatatip(n),'DataCursor'),'DataIndex',index(n, :),...
      'TargetPoint',pos)
  end
  
  % Specify datatip properties
  set(hDatatip(n),'Position',pos)
  
end

% Update all data cursors
updateDataCursors(hDataCursorMgr)

% Return handles if requested
if nargout==1
  varargout = {hDatatip};
end

end
