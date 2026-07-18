% fsle - calculate the finite size lyapunov exponent
% lamba = 1/tau*log(df/d0);

%% INITIAL SETTINGS
d0 = delta; % initial, uniform displacement from primary to secondary seeds
df = 3*d0; % threshold distance, this is what makes FSLE finite-SIZE
C = log(df/d0)

disp('FSLE - calculate finite size Lyapunov exponent.');

%% FIND TIME STEP AT WHICH DF IS REACHED FOR EACH SECONDARY SEED
for itdf = 1:size(distUp,1);
dfUp = find(distUp(itdf,:) >= df);
if isempty(dfUp);
dfUpFirst(itdf,1) = NaN;
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
dfDnFirst(itdf,1) = NaN;
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
dfLtFirst(itdf,1) = NaN;
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
dfRtFirst(itdf,1) = NaN;
else if ~isempty(dfRt);
dfRtFirst(itdf,1) = dfRt(1);
dfRtFirstVal(itdf,1) = distRt(itdf,dfRt(1));
    end
end
clear dfRt
end

% combine all df datasets into one searchable data set
dfAll = [dfUpFirst,dfDnFirst,dfLtFirst,dfRtFirst];
minDfAll = min(dfAll,[],2);
tauDfAll = minDfAll.*dt;
%%


%% CALCULATE FSLE
lambdaFSLE = (1./tauDfAll).*log(df/d0);
%%

%% CALCULATE FTLE?
%T = 500; % time step to calculate at
%dfT = [distUp(:,T),distDn(:,T),distLt(:,T),distRt(:,T)];
%dfT = max(dfT,[],2);
%lambdaFTLE = 1/T * log(dfT./d0);
%%

%% PLOTS
% Plot sorted max displacements for each secondary seed
%figure();
%plot(sort(maxUp),'r');
%hold on; %axis equal;
%plot(sort(maxDn),'g');
%plot(sort(maxLt),'c');
%plot(sort(maxRt),'b');

figure();
scatter(xn(:,1),yn(:,1),8,lambdaFSLE,'filled');
view(2);
axis equal
colorbar

%%
disp('FSLE field calculated.');