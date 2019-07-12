% flowmap2d - create 2D (xy) flow map from particle tracers
disp('FLOWMAP2D - visualize flow map in xy plane.');

close all

% requires output from ADVECT
contourFlowmap = 1;

%% 
% calculate displacements between time steps
dxnS = xnS(:,2:end) - xnS(:,1:end-1);
dynS = ynS(:,2:end) - ynS(:,1:end-1);

% calculate euclidean 2d distance between time steps
distXYn = bsxfun(@hypot,dxnS,dynS);
for itFlowMap = 1:TS;
    flowmap(:,itFlowMap) = sum(distXYn(:,1:itFlowMap),2);
end % end loop over all TS
%%

%% INTERPOLATE TO GRID
uFlowmap = TriScatteredInterp(xnS(:,1),ynS(:,2),flowmap(:,TS),'natural');
gridFlowmap = uFlowmap(seedXBack,seedYBack);
%%

%% PLOT FLOW MAP
%n = 100;
%scatter(xnS(:,1),ynS(:,1),10,flowmap(:,TS),'filled');
%colorbar;
%axis equal;
%%

if contourFlowmap == 1;
% CONTOUR PLOT VELOCITY FIELD
    figure();
    contourf(seedXBack,seedYBack,gridFlowmap,10);
    colorbar
    axis equal
    hold on;
    plot(seedXY(:,1),seedXY(:,2),'.k');
    axis equal
    xlabel('x [m]');
    ylabel('y [m]');
    legend('\phi [m]','initial seed locations');
else
end

figure();
plot(distXYn(6,:));
