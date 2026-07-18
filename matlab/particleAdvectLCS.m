% particleAdvect
% need to add nearest interpolant, current version returns too many points!
% STEPS:
% 1) Export nodes and velocity field data into text files, corresponding to
% nodal locations
% 2) Create a velocity.mat file for the velocity field data
% 3) Import the nodal xyz locations into superScat
% 4) Run particleAdvect
% 5) Run superScat to export advection data set to SMS

% Goals: increase advection speed for high-res sets

% Note: >"res" affects speed by increasing the total number of seed locations
%       >the size of the initial velocity field affects speed by requiring
%       more locations to be searched and interpolated

%clear all
%close all

%% Import velocity and node scatter data
%load('F:\Tuhtan Laboratory\MATLAB Test\superscatter\rottenburg\velocity.mat');
%load('F:\Tuhtan Laboratory\MATLAB Test\superscatter\rottenburg\xyz.mat');
%load('H:\superscatter\rottenburg\velocity.mat');
%load('H:\superscatter\rottenburg\xyz.mat');
%

% load presaved gyrefield 
load gyrefieldOUT.mat
disp('Velocity and nodal scatter data loaded.');

%% Initial settings
defaultElev = 0; % default elevation
numScatter = size(xyz,1); % total number of node scatter points
TS = 50; % number of time steps to include
res = 0.01; % high res fish pass 0.10, Rottenburg 3.00
visualSMS = 0; % prepare 'out' matrix for SMS visualization
%%

dt = 10.00; % high res fish pass dt = 0.025
n = TS;
xOffset = 8;
delta0 = 1*gridRes; % (2*gridRes a good default) offset distance from initial particles at [x0, y0]

xMin = min(xyz(:,1))+2*delta0; xMax = 1.00*max(xyz(:,1))-2*delta0; % explicit min, max is not used because offset NN particles will be excluded!
yMin = min(xyz(:,2))+2*delta0; yMax = 1.00*max(xyz(:,2)-2*delta0); % 2*delta0 seems to work well in SMS visualization

[SeedX0, SeedY0] = meshgrid(xMin:res:xMax,yMin:res:yMax);

groupsOn = 1;

if groupsOn == 1;

sizeGroup = size(SeedX0,1);    
    
% Create displaced seeds 
SeedXA = SeedX0 + delta0;
SeedYA = SeedY0;

SeedXB = SeedX0 - delta0;
SeedYB = SeedY0;

SeedXC = SeedX0;
SeedYC = SeedY0 + delta0;

SeedXD = SeedX0;
SeedYD = SeedY0 - delta0;

SeedX = [SeedXA;SeedXB;SeedXC;SeedXD];
SeedY = [SeedYA;SeedYB;SeedYC;SeedYD];

% Assign constant values to each particle group
val0 = 1;
valA = 2;
valB = 3;
valC = 4;
valD = 5;
else
end

disp('Seed particles initialized.');

% Create interpolation surfaces based on x y u and x y v
% method can be either natural, linear, nearest
% natural has previously caused MATLAB to crash
uF = TriScatteredInterp(xyz(:,1),xyz(:,2),velocity(:,1),'nearest'); % velocity and x y datasets must correspond to the same nodes!
vF = TriScatteredInterp(xyz(:,1),xyz(:,2),velocity(:,2),'nearest');
disp('Velocity field surfaces complete.');

% Evaluate velocities at all grid locations
Seedu = uF(SeedX,SeedY);
Seedv = vF(SeedX,SeedY);

% Evaluate velocities at initial positions only
Seedu0 = uF(SeedX0,SeedY0);
Seedv0 = vF(SeedX0,SeedY0);

% Assign 0 values to NaN
idxuSeed = find(Seedu==0);
Seedu(idxuSeed) = NaN;

idxuSeed0 = find(Seedu0==0);
Seedu0(idxuSeed0) = NaN;

% Calculate velocity magnitude at each seed location
magSeed = (Seedu.^2 + Seedv.^2).^0.5;
magSeed0 = (Seedu0.^2 + Seedv0.^2).^0.5;
disp('Seed velocity magnitudes calculated.');

%figure;
%surf(magSeed);
%axis equal
%shading flat

idxNaN = find(~isnan(magSeed));
un(:,1) = Seedu(idxNaN);
vn(:,1) = Seedv(idxNaN);
mag(:,1) = (un(:,1).^2 + vn(:,1).^2).^0.5;
idxMagMin = find(mag(:,1)<=0.01);
mag(idxMagMin,1) = 0.01;

idx0NaN = find(~isnan(magSeed0));
u0n(:,1) = Seedu0(idx0NaN);
v0n(:,1) = Seedv0(idx0NaN);
mag0(:,1) = (u0n(:,1).^2 + v0n(:,1).^2).^0.5;
idxMag0Min = find(mag0(:,1)<=0.01);
mag0(idxMag0Min,1) = 0.01;

if groupsOn == 1;
group = zeros(size(SeedX,1),size(SeedX,2)); % initialize group array
group(1:sizeGroup,:) = valA;
group(sizeGroup+1:2*sizeGroup,:) = valB;
group(2*sizeGroup+1:3*sizeGroup,:) = valC;
group(3*sizeGroup+1:4*sizeGroup,:) = valD;
groupOut(:,1) = group(idxNaN); 
else
end

xn(:,1) = SeedX(idxNaN);
yn(:,1) = SeedY(idxNaN);

x0n(:,1) = SeedX0(idx0NaN);
y0n(:,1) = SeedY0(idx0NaN);

numSeeds = size(xn,1);
num0 = size(x0n,1);

% advect seed particle locations
disp('Begin seed particle advection.')
for it = 2:n;
    % advect all positions
    xn(:,it) = xn(:,it-1) + dt*un(:,it-1);
    yn(:,it) = yn(:,it-1) + dt*vn(:,it-1);
    disp(['Group positions advected ', num2str(it) '/', num2str(n)]);
    % evaluate new velocities
    un(:,it) = uF(xn(:,it-1),yn(:,it-1));
    vn(:,it) = vF(xn(:,it-1),yn(:,it-1));
    disp(['Group velocities updated ', num2str(it) '/', num2str(n)]);
    mag(:,it) = (un(:,it).^2 + vn(:,it).^2).^0.5;
    idxMagMin = find(mag(:,it)<=0.01); % to keep pts from vanishing in SMS, keep some minimum velmag value for all "slow" pts
    mag(idxMagMin,it) = 0.01;
    
    % advect only initial positions x0, y0
    x0n(:,it) = x0n(:,it-1) + dt*u0n(:,it-1);
    y0n(:,it) = y0n(:,it-1) + dt*v0n(:,it-1);
    disp(['Initial seed positions advected ', num2str(it) '/', num2str(n)]);
    % evaluate new velocities
    u0n(:,it) = uF(x0n(:,it-1),y0n(:,it-1));
    v0n(:,it) = vF(x0n(:,it-1),y0n(:,it-1));       
    disp(['Initial seed velocities updated ', num2str(it) '/', num2str(n)]);
    %plot(xn(:,it),yn(:,it),'.b');
end % end loop over all time steps

idxMagNaN = find(isnan(mag));
mag(idxMagNaN) = 0;
xn(idxMagNaN) = 0;
yn(idxMagNaN) = 0;

%idxSelectA = find(xyz(:,1)<486438); % find all initial locations less than fixed x value
%idxSelectB = find(xyz(:,1)>9); % find all initial locations greater than fixed x value

%% CREATE DATA SET OUT FOR SMS VISUALIZATION
% This routine requires the mapping (interpolation) of the actual particle locations to the
% corresponding fixed positions of the super scatter data set for subsequent visualization in SMS.
% This loop is NOT necessary if the results are not to be visualized in SMS!
if visualSMS == 1;
out = zeros(numScatter,TS); % initialze matrix out
for it = 1:n;
% create dataset xnyn for each time step, single array with x y coords of advected points
xnyn(:,1) = xn(:,it);
xnyn(:,2) = yn(:,it);

% index nearest corresponding points on the fixed scatter map xyz
% xnyn are the "real" positions of the advected particles
% xyz(:,1:2) are the "fixed" positions of the underlying nodal scatter set
% NN creates a mapping from the real to the fixed positions
NN = ipdm(xnyn,xyz(:,1:2),'Subset','nearest','Result','Structure'); % NOTE: ipdm.m is a dependency!
idxNN = NN.columnindex;
%out(idxNN,it) = 10; % this sets the default value of all visualized super scatter points
%out(idxNN',it) = mag(:,it); % assign pt elevation as vel mag
out(idxNN',it) = groupOut(:,1); % assign pt elevation as constant value for each group

%idxTest = find(xn(:,1)<6);
%NNselect = ipdm(xnyn(idxTest,:),xyz(:,1:2),'Subset','nearest'); % NOTE: ipdm.m is a dependency!
%idxNNselect = find(any(NNselect)); % index nearest points in xyz
%out(idxNNselect,it) = 20;

%membersA = ismember(idxSelectA,idxNN); % finds common ids in both sets
%idxA = idxSelectA(membersA);
%out(idxA,it) = 10;
%clear membersA

%membersB = ismember(idxSelectB,idxNN);
%idxB = idxSelectB(membersB);
%out(idxB,it) = 2;
%clear membersB

display(['TS ' num2str(it) '/' num2str(n) ' complete.']);

end % end loop over all time steps n
else
end % do not run this loop if SMS visualization is not required
%%


%% PLOTS
figure;
plot(xn,yn,'.');
axis equal
hold on
plot(x0n,y0n,'.k');
%%