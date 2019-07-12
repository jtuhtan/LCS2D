% flowmap - create a 2D flow map from a 2D velocity field
disp('FLOWMAP - create a flow map from a 2D velocity field');

%% 1 USE "SEED" TO SEED FLOW FIELD
%seed
%%

%% 2 ADVECT SEEDS
%advect
%%

%% 3 CALCULATE EUCLIDEAN DISTANCE FROM X0 TO XTn
 ynDelta = bsxfun(@minus,yn,yn(:,1));
 xnDelta = bsxfun(@minus,xn,xn(:,1));
 
 eucDist = sqrt(xnDelta.^2 + ynDelta.^2);
 %% 
 
 %% PLOTS
 % plot of all pts euclidean distance as a function of time
 plot(eucDist,1:TS+1,'.')
 
 %%