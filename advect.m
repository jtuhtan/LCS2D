% advect - advects particle seeds from a preloaded flow map and primary seeds
disp('ADVECT - advect particle seeds.');

%%
% Jeff Tuhtan 2013 - Distributable GPL
% jtuhtan@gmail.com
% V.20130208
%%

%%
% DEPENDENCIES: SEED 
% TODO: find interpolant faster than MATLAB's TriScatteredInterp
% evaluation? - causes major bottleneck in advection routine.
% TODO: add interpoint distance matrix to check if max eig alignment occurs
% for small delta. use alpha_i = delta_n / delta_0 where i is the secondary
% seed. If alpha is more or less the same for all seeds, then only one
% secondary seed need be advected to approx. FTLE!fi
%%

%%
% REQUIRED INPUTS: seedXY, seedUV, resSeed
%%

close all

%% INITAL SETTINGS
plotsOn = 1;
TS = 10; % number of internal time steps for a stationary flow map
dt = 0.02; % size of advection time step
delta = 0.001*resSeed; % fixed secondary perturbation distance x y from the primary seed
% NOTE: smaller dt means that the advection will be truer to the flow map,
% but will require larger values of TS for the same total advection time.
% In order to keep the mapping accurate, it may be best to decrease the
% seed resolution rather than to increase dt or decrease TS.
timeAdvect = TS * dt; % total time particles are advected for a stationary flow map
disp(['No. time steps, TS: ', num2str(TS), ' sec']);
disp(['Advection time step, dt: ', num2str(dt), ' sec']);
disp(['Total advection time for stationary flow map: ', num2str(timeAdvect), ' sec']);
%%

%% LOAD PARTICLE SEEDS
% load presaved seed locations at initial positions (x0,y0,t0)
numSeed = size(seedXY,1); % get total no. of seeds
disp('Initial seed positions and velocities loaded.');
%%

%% CHECK CFL CONDITION IN X AND Y
uMax = max(seedUV(:,1));
vMax = max(seedUV(:,2));
velMax = max(uMax,vMax);
seedRes = ((seedXY(2,1) - seedXY(1,1))^2 + (seedXY(2,2) - seedXY(2,1))^2)^0.5;
velAdvect = seedRes / dt; 
CFLMax = 1; % set threshold for CFL criteria
CFL = uMax / velAdvect + vMax / velAdvect;

disp(['Initial CFL [<1]: ', num2str(CFL)]);
if CFL > 1;
    disp('CFL criteria not met! Recommend reducing seed resolution or dt. Aborting advection.');
else
    disp('CFL criteria satisfied.');
end
    
if CFL <= 1;
%%

%% ADVECT PRIMARY SEEDS
xnS = seedX(:);
ynS = seedY(:);
% intialize velocities
unS(:,1) = uF(xnS(:,1),ynS(:,1));
vnS(:,1) = vF(xnS(:,1),ynS(:,1));

for itTS = 2:TS+1;
    xnS(:,itTS) = xnS(:,itTS-1) + dt*unS(:,itTS-1);
    ynS(:,itTS) = ynS(:,itTS-1) + dt*vnS(:,itTS-1);
    % get new velocities
    unS(:,itTS) = uF(xnS(:,itTS-1),ynS(:,itTS-1));
    vnS(:,itTS) = vF(xnS(:,itTS-1),ynS(:,itTS-1));
    disp(['Advected Primary Seeds: ', num2str(itTS-1) '/' num2str(TS)]);
end % end loop over all TS
%%

%% TEST SINGLE SECONDARY
testSingleSecondary = 0;
if testSingleSecondary == 1;
xnS = zeros(numSeed,TS);
ynS = xnS;
xnS(:,1) = xn(:,1) + 0; 
ynS(:,1) = yn(:,1) + delta;
% intialize velocities
unS(:,1) = uF(xnS(:,1),ynS(:,1));
vnS(:,1) = vF(xnS(:,1),ynS(:,1));

    for itTS = 2:TS+1;
    xnS(:,itTS) = xnS(:,itTS-1) + dt*unS(:,itTS-1);
    ynS(:,itTS) = ynS(:,itTS-1) + dt*vnS(:,itTS-1);
    % get new velocities
    unS(:,itTS) = uF(xnS(:,itTS-1),ynS(:,itTS-1));
    vnS(:,itTS) = vF(xnS(:,itTS-1),ynS(:,itTS-1));
    disp(['Advected Secondary S: ', num2str(itTS-1) '/' num2str(TS)]);
    
    % Plot primary advection
    %hold on;
    %plot(xnS(:,itTS),ynS(:,itTS),'.r');
    %
    end % end loop over all TS
else
end
%%

%% ADVECT: SECONDARY SEEDS EULER ONE-STEP
%             UP
%           delta
% LT delta (xnS,ynS) delta RT
%           delta
%             DN
advectSecondary = 1;
if advectSecondary == 1;
% initialize matrices
xnUp = zeros(numSeed,TS);
ynUp = xnUp;

xnDn = xnUp;
ynDn = xnUp;

xnLt = xnUp;
ynLt = xnUp;

xnRt = xnUp;
ynRt = xnUp;

% offset from primary particles at t0
xnUp(:,1) = xnS(:,1) + 0; ynUp(:,1) = ynS(:,1) + delta;
xnDn(:,1) = xnS(:,1) + 0; ynDn(:,1) = ynS(:,1) - delta;
xnLt(:,1) = xnS(:,1) - delta; ynLt(:,1) = ynS(:,1) + 0;
xnSRt(:,1) = xnS(:,1) + delta; ynRt(:,1) = ynS(:,1) + 0;
%

Up = [0 delta];
Dn = [0 -delta];
Lt = [-delta 0];
Rt = [delta 0];

OffsetNames = {'Up';'Down';'Left';'Right'};

Offsets = [Up; Dn; Lt; Rt]; % combines all offsets into one matrix

xnOff = zeros(numSeed,itTS);
ynOff = zeros(numSeed,itTS);
unOff = zeros(numSeed,itTS);
vnOff = zeros(numSeed,itTS);

for itOff = 1:4;
    % create initial matrix using offsets
    xnOff(:,1) = xnS(:,1)+Offsets(itOff,1); 
    ynOff(:,1) = ynS(:,1)+Offsets(itOff,2);
    unOff(:,1) = uF(xnOff(:,1),ynOff(:,1));
    vnOff(:,1) = vF(xnOff(:,1),ynOff(:,1));
for itTS = 2:TS+1;
    xnOff(:,itTS) = xnOff(:,itTS-1) + dt*unOff(:,itTS-1);
    ynOff(:,itTS) = ynOff(:,itTS-1) + dt*vnOff(:,itTS-1);
    % get new velocities
    unOff(:,itTS) = uF(xnOff(:,itTS),ynOff(:,itTS));
    vnOff(:,itTS) = vF(xnOff(:,itTS),ynOff(:,itTS));
    
    disp(['Advected Secondary ' , OffsetNames{itOff,1}, ': ' num2str(itTS-1) '/' num2str(TS)]);
end % end loop over all TS

if itOff == 1;
    xnUp = xnOff; ynUp = ynOff; unUp = unOff; vnUp = vnOff;
else if itOff == 2;
        xnDn = xnOff; ynDn = ynOff; unDn = unOff; vnDn = vnOff;
    else if itOff == 3;
            xnLt = xnOff; ynLt = ynOff; unLt = unOff; vnLt = vnOff;
        else if itOff == 4;
                xnRt = xnOff; ynRt = ynOff; unRt = unOff; vnRt = vnOff;
            else
            end
        end
    end
end  
    
end % end loop over all offsets
else
end
%%

%% PLOTS
if plotsOn == 1;
% Plot initial seed locations, primary and secondary
figure();
plot(xnS(:,1),ynS(:,1),'+k');
hold on; axis equal;
plot(xnUp(:,1),ynUp(:,1),'.r');
plot(xnDn(:,1),ynDn(:,1),'.g');
plot(xnLt(:,1),ynLt(:,1),'.c');
plot(xnRt(:,1),ynRt(:,1),'.b');

% Plot all seed locations
figure();
plot(xnS,ynS,'.k');
hold on; axis equal;
plot(xnUp,ynUp,'.r');
plot(xnDn,ynDn,'.g');
plot(xnLt,ynLt,'.c');
plot(xnRt,ynRt,'.b');

% Plot only final seed locations
figure();
plot(xnS(:,TS),ynS(:,TS),'+k');
hold on; axis equal;
plot(xnUp(:,TS),ynUp(:,TS),'.r');
plot(xnDn(:,TS),ynDn(:,TS),'.g');
plot(xnLt(:,TS),ynLt(:,TS),'.c');
plot(xnRt(:,TS),ynRt(:,TS),'.b');
else
end
%%
disp('All particles advected.');
disp(['Total advection time: ', num2str(dt*TS), ' seconds']);
else
end