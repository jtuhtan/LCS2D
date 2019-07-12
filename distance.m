% distance - calculates distances between primary and secondary seeds
disp('DISTANCE - calculate distances between all primary and secondary particle seeds.');

%% CALCULATE STRETCHING -"UP"
deltaXUp = xn - xnUp; % calculate pointwise difference between xn and xnUp
deltaYUp = yn - ynUp; % calculate pointwise difference between yn and ynUp
distUp = (deltaXUp.^2 + deltaYUp.^2).^0.50; % calculate euclidean distance 2D between primary and secondary "up"
%%

%% CALCULATE STRETCHING -"DOWN"
deltaXDn = xn - xnDn; % calculate pointwise difference between xn and xnUp
deltaYDn = yn - ynDn; % calculate pointwise difference between yn and ynUp
distDn = (deltaXDn.^2 + deltaYDn.^2).^0.50; % calculate euclidean distance 2D between primary and secondary "down"
%%

%% CALCULATE STRETCHING -"LEFT"
deltaXLt = xn - xnLt; % calculate pointwise difference between xn and xnUp
deltaYLt = yn - ynLt; % calculate pointwise difference between yn and ynUp
distLt = (deltaXLt.^2 + deltaYLt.^2).^0.50; % calculate euclidean distance 2D between primary and secondary "left"
%%

%% CALCULATE STRETCHING -"RIGHT"
deltaXRt = xn - xnRt; % calculate pointwise difference between xn and xnUp
deltaYRt = yn - ynRt; % calculate pointwise difference between yn and ynUp
distRt = (deltaXRt.^2 + deltaYRt.^2).^0.50; % calculate euclidean distance 2D between primary and secondary "right"
%%
disp('All distances calculated.');

