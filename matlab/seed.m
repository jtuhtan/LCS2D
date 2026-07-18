% seed - creates initial particle seeds for a preloaded flow map
disp('SEED - create particle seeds.');

%%
% Jeff Tuhtan 2013 - Distributable GPL
% jtuhtan@gmail.com
% V.20130208
%%

clear all
close all

%% LOAD PRESAVED SCATTER AND FLOW MAP
%load gyrefieldOUT.mat
load GertData.mat
% need velocity [u v], xyz [x y z]
disp('Velocity and nodal scatter data loaded.');
%%

%% SET SEED GRID RESOLUTION
resSeed = 10; % resolution of seed spacing
interpMethod = 'natural'; % interpolation method used for velocity field (natural, linear, nearest)
resBack = resSeed;%resGridVel; % resolution of background contour plot, only used for visualization
%%

%% CALCULATE VELOCITY MAGNITUDE
velocityMag = bsxfun(@hypot,velocity(:,1),velocity(:,2));
%%

%% PLOT SCATTER NODES & VELOCITY DATA
figure();
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),3,xyz(:,3),'filled');
axis equal;
figure();
scatter(xyz(:,1),xyz(:,2),3,velocityMag,'filled');
colorbar
axis equal;
%%

%% INITIAL SETTINGS
magSeedThresh = 0; % threshold minimum velocity magnitude used to remove initial seeds
disp(['Uniform seed spacing: ', num2str(resSeed), ' m']);
disp(['Minimum threshold velocity: ', num2str(magSeedThresh), ' m/s']);
%%

%% PERFORM INITIAL CALCULATIONS
% retrieve scatter extents in x y plane
offset = 2*resSeed; % offset distance from absolute boundary of the inital dataset, helps remove spurious interpolation at borders

% retrieve absolute limits of xy min, max from scatter set
xMinAbs = min(xyz(:,1)); xMaxAbs = max(xyz(:,1));
yMinAbs = min(xyz(:,2)); yMaxAbs = max(xyz(:,2));

% create seed offset xy min, max
xMin = min(xyz(:,1)) + offset; xMax = max(xyz(:,1)) - offset;
yMin = min(xyz(:,2)) + offset; yMax = max(xyz(:,2)) - offset;
disp(['xMin: ', num2str(xMinAbs), ' xMax: ', num2str(xMaxAbs)]);
disp(['yMin: ', num2str(yMinAbs), ' yMax: ', num2str(yMaxAbs)]);

% OR - manually define min max values
%xMin = 7; xMax = 69;
%yMin = 5; yMax = 31; 

% create uniformly distributed seeds covering entire extents in x y plane
[seedX, seedY] = meshgrid(xMin:resSeed:xMax,yMin:resSeed:yMax);
[seedXBack, seedYBack] = meshgrid(xMin:resBack:xMax,yMin:resBack:yMax);
xMinSeed = min(seedX(:)); xMaxSeed = max(seedX(:));
yMinSeed = min(seedY(:)); yMaxSeed = max(seedY(:));
xyArea = (xMaxSeed - xMinSeed) * (yMaxSeed - yMinSeed);

xyAreaTotal = (xMaxAbs - xMinAbs) * (yMaxAbs - yMinAbs);
xySeedDensity = (size(seedX,1) * size(seedX,2)) / xyArea;
xyScatterDensity = (size(xyz,1)) / xyArea;
xyDensityRatio = xySeedDensity / xyScatterDensity;
xyMinRes = resSeed * sqrt(xyDensityRatio);

disp(['xy planar area: ', num2str(xyAreaTotal), ' m^2']);
disp(['Scatter density: ', num2str(ceil(xyScatterDensity)), ' pts/m^2']);
disp(['Seed density: ', num2str(ceil(xySeedDensity)), ' pts/m^2']);
disp(['Ratio Seed/Scatter [should be > 1]: ', num2str(xyDensityRatio)]);
disp(['Minimum seed resolution: ', num2str(xyMinRes)]);
% Create interpolation surfaces based on x y u and x y v
% method can be either natural, linear, nearest
% NOTE: natural has previously caused MATLAB to crash!
disp('Interpolating velocity fields...');
uF = TriScatteredInterp(xyz(:,1),xyz(:,2),velocity(:,1),interpMethod); % velocity and x y datasets must correspond to the same nodes!
vF = TriScatteredInterp(xyz(:,1),xyz(:,2),velocity(:,2),interpMethod);
disp('Velocity field surfaces complete.');

% Evaluate velocities at all grid locations
seedU = uF(seedX,seedY);
seedV = vF(seedX,seedY);

seedUBack = uF(seedXBack,seedYBack);
seedVBack = vF(seedXBack,seedYBack);

% Calculate velocity magnitude at each seed location
magSeed = sqrt(seedU.^2 + seedV.^2);
magSeedBack = sqrt(seedUBack.^2 + seedVBack.^2);

idxMagSeed = find(magSeed <= magSeedThresh);
idxMagSeedKeep = find(magSeed > magSeedThresh);
disp('Seed velocity magnitudes calculated.');
%%

%% FIX SPURIOUS VALUES
velocityMag(isnan(velocityMag)) = 0;
badValue = 0; % manual assignment of spurious value
idxBad = find(velocityMag==badValue); % index bad values
idxGood = find(~velocityMag==badValue);
dataGood = xyz(idxGood,:); % create filtered data set containing only good values
velocityGood = velocity(idxGood,:);

uFGood = TriScatteredInterp(dataGood(:,1),dataGood(:,2),velocityGood(:,1),'natural'); % velocity and x y datasets must correspond to the same nodes!
vFGood = TriScatteredInterp(dataGood(:,1),dataGood(:,2),velocityGood(:,2),'natural');

uF = uFGood;
vF = vFGood;

% Evaluate velocities at all grid locations
seedUGood = uFGood(seedX,seedY);
seedVGood = vFGood(seedX,seedY);

% Calculate velocity magnitude at each seed location 
magSeed = sqrt(seedUGood.^2 + seedVGood.^2);
idxMagSeed = find(magSeed <= magSeedThresh);
idxMagSeedKeep = find(magSeed > magSeedThresh);
disp('Seed velocity magnitudes calculated.');
%%

%% OPTIONAL CALCULATIONS
replyMagFilter = 0;
disp(['Number of seeds at or below minimum velocity threshold: ', num2str(size(idxMagSeed,1))]);
reply = input('Filter seeds using minimum velocity threshold? Y/N [N]: ', 's' );
if strcmp(reply,'y');
    replyMagFilter = 1;
        else if strcmp(reply,'Y');
                    replyMagFilter = 1;
                else
                    disp('Seed velocities NOT filtered.');
            end
end
% assign the filtered data to a constant value
if replyMagFilter == 1;
    seedX(idxMagSeed) = NaN;
    seedY(idxMagSeed) = NaN;
    magSeed(idxMagSeed) = NaN;
    disp('Seed velocities filtered.');
else
end 
%%

%% CREATE COLUMN VECTOR WITH ALL XY SEED LOCATIONS
seedXY = [seedX(:),seedY(:)];
seedUV = [seedUGood(:),seedVGood(:)];
disp('X Y seeds created.');
%%

%% PLOT SEEDS
replyPlotY = 0;
replyPlot = input('Plot seeds? Y/N [N]: ', 's');
if strcmp(replyPlot,'y');
    replyPlotY = 1;
else if strcmp(replyPlot,'Y');
        replyPlotY = 1;
    end
end

if replyPlotY == 1;
% CONTOUR PLOT VELOCITY FIELD
    figure();
    contourf(seedXBack,seedYBack,magSeedBack,10);
    colorbar
    axis equal
    hold on;
    plot(seedXY(:,1),seedXY(:,2),'.k');
    axis equal
    xlabel('x [m]');
    ylabel('y [m]');
    legend('velocity magnitude [m/s]','primary seed locations');
else
end
%%

