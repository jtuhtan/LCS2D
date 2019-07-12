% fsle - calculate the finite size lyapunov exponent
% lamba = 1/tau*log(df/d0);

%% INITIAL SETTINGS
d0 = delta; % initial, uniform displacement from primary to secondary seeds
df = 2*d0; % threshold distance, this is what makes FSLE finite-SIZE
%%

disp('FSLE - calculate finite size Lyapunov exponent.');

%% FIND TIME STEP AT WHICH DF IS REACHED FOR EACH SECONDARY SEED
for itdf = 1:size(distUp,1);
dfUp = find(distUp(itdf,:) >= df);
if isempty(dfUp);
dfUpFirst(itdf,1) = 0;
else if ~isempty(dfUp);
dfUpFirst(itdf,1) = dfUp(1);
dfUpFirstVal(itdf,1) = distUp(itdf,dfUp(1));
    end
end
clear dfUp
end

for itdf = 1:size(distUp,1);
dfDn = find(distDn(itdf,:) >= df);
if isempty(dfDn);
dfDnFirst(itdf,1) = 0;
else if ~isempty(dfDn);
dfDnFirst(itdf,1) = dfDn(1);
dfDnFirstVal(itdf,1) = distDn(itdf,dfDn(1));
    end
end
clear dfDn
end

for itdf = 1:size(distUp,1);
dfLt = find(distLt(itdf,:) >= df);
if isempty(dfLt);
dfLtFirst(itdf,1) = 0;
else if ~isempty(dfLt);
dfLtFirst(itdf,1) = dfLt(1);
dfLtFirstVal(itdf,1) = distLt(itdf,dfLt(1));
    end
end
clear dfLt
end

for itdf = 1:size(distUp,1);
dfRt = find(distRt(itdf,:) >= df);
if isempty(dfRt);
dfRtFirst(itdf,1) = 0;
else if ~isempty(dfRt);
dfRtFirst(itdf,1) = dfRt(1);
dfRtFirstVal(itdf,1) = distRt(itdf,dfRt(1));
    end
end
clear dfRt
end

% combine all df datasets into one searchable data set
dfAll = [dfUpFirst,dfDnFirst,dfLtFirst,dfRtFirst];
tau

%%

%% FIND MAXIMUM DISTANCES FOR EACH SECONDARY SEED
maxUp = max(distUp,[],2);
maxDn = max(distDn,[],2);
maxLt = max(distLt,[],2);
maxRt = max(distRt,[],2);

for itMax = 1:size(distUp,1);
    idxMaxVal = find(distUp(itMax,:)==maxUp(itMax,1)); % returns FIRST max value, in case the same value occurs more than once
    if size(idxMaxVal,2) > 1;
        idxMaxVal = idxMaxVal(1);
    else
    end
    idxMaxFirst(itMax,1) = idxMaxVal;
    tauUp(itMax,1) = idxMaxVal*dt;
    tsUp(itMax,1) = idxMaxVal;
end % end loop over all seed points

for itMax = 1:size(distDn,1);
    idxMaxVal = find(distDn(itMax,:)==maxDn(itMax,1)); % returns FIRST max value, in case the same value occurs more than once
    if size(idxMaxVal,2) > 1;
        idxMaxVal = idxMaxVal(1);
    else
    end
    idxMaxFirst(itMax,1) = idxMaxVal;
    tauDn(itMax,1) = idxMaxVal*dt;
    tsDn(itMax,1) = idxMaxVal;
end % end loop over all seed points

for itMax = 1:size(distLt,1);
    idxMaxVal = find(distLt(itMax,:)==maxLt(itMax,1)); % returns FIRST max value, in case the same value occurs more than once
    if size(idxMaxVal,2) > 1;
        idxMaxVal = idxMaxVal(1);
    else
    end
    idxMaxFirst(itMax,1) = idxMaxVal;
    tauLt(itMax,1) = idxMaxVal*dt;
    tsLt(itMax,1) = idxMaxVal;
end % end loop over all seed points

for itMax = 1:size(distRt,1);
    idxMaxVal = find(distRt(itMax,:)==maxRt(itMax,1)); % returns FIRST max value, in case the same value occurs more than once
    if size(idxMaxVal,2) > 1;
        idxMaxVal = idxMaxVal(1);
    else
    end
    idxMaxFirst(itMax,1) = idxMaxVal;
    tauRt(itMax,1) = idxMaxVal*dt;
    tsRt(itMax,1) = idxMaxVal;
end % end loop over all seed points
%%

%% CALCULATE TAU (TIME TO MAX DISPLACEMENT) FOR EACH SECONDARY SEED
% we need to select the maximum displacement, and its corresponding time
tsAll = [tsDn,tsUp,tsLt,tsRt];
maxAll = [maxUp,maxDn,maxLt,maxRt];
tauAll = [tauUp,tauDn,tauLt,tauRt];
maxAllVal = max(maxAll,[],2);
for itMax = 1:size(distRt,1);
    getMaxAllVal = find(maxAll(itMax,:)==maxAllVal(itMax,1));
    if size(getMaxAllVal,2)>1;
        getMaxAllVal = getMaxAllVal(1); % retrieves first max value
    else
    end
    idxMaxAllTS(itMax,1) = getMaxAllVal; 
    tsMax(itMax,1) = tsAll(idxMaxAllTS(itMax,1));
end
tauMax = tsMax.*dt;

%%

%% CALCULATE FSLE
lambdaFSLE = (1./tauMax).*log(maxAllVal./delta);
%%

%% PLOTS
% Plot sorted max displacements for each secondary seed
figure();
plot(sort(maxUp),'r');
hold on; %axis equal;
plot(sort(maxDn),'g');
plot(sort(maxLt),'c');
plot(sort(maxRt),'b');

figure();
scatter(xn(:,1),yn(:,1),8,lambdaFSLE,'filled');
axis equal
colorbar

%%
disp('FSLE field calculated.');