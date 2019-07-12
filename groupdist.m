% group distance - calculates the distances between the initial particle
% and its corresponding group

% Dependencies - particleAdvectLCS, ipdm

% Jeff A. Tuhtan 2013
% jtuhtan@gmail.com

%% 1 - Locate the nearest n neigbors based on initial (t0,x0,y0) positions

for itGroupDist = 1:size(x0n,1);
    
SeedsInitial = [xn(:,1),yn(:,1)]; % grabs all seeds for initial time step 
PointsInitial = [x0n(itGroupDist,1),y0n(itGroupDist,1)]; % grabs one coordinate pair at initial time step

NNgroup = ipdm(PointsInitial,SeedsInitial,'Subset','smallestfew','limit',4,'Result','Structure'); % NOTE: ipdm.m is a dependency!
idxGroupNN(:,itGroupDist) = NNgroup.columnindex;

disp([num2str(itGroupDist) ' / ' num2str(size(x0n,1))]);
end % end loop over all initial points

%%

%% 2 - Find distances between all points and for all time steps

for itGroupTS = 1:TS;
for itGroupDist = 1:size(x0n,1);
    
SeedsTS = [xn(idxGroupNN(:,itGroupDist),itGroupTS),yn(idxGroupNN(:,itGroupDist),itGroupTS)]; % Grabs the new locations of the seeds
PointsTS = [x0n(itGroupDist,itGroupTS),y0n(itGroupDist,itGroupTS)]; % Grabs a single point    
    
distGroup(itGroupDist,:,itGroupTS) = ipdm(PointsTS,SeedsTS); % returns the distances between the original nearest neighbors over each time step

distMaxGroup(itGroupDist,itGroupTS) = max(distGroup(itGroupDist,:,itGroupTS)); % returns only the maximum distance of the original nearest neighbors
end % end loop over all initial points
disp([num2str(itGroupTS) ' / ' num2str(TS)]);
end % end loop over all TS

%%