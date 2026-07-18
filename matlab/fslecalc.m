%% fslecalc - calculates the finite size Lyapunov exponent 
% lamda = 1/tau * log(df/d0)
% tau - time required to reach df
% d0 - initial seed neighbor displacement
% df - threshold seed neighbor displacement
% NOTE: if d0 and df are a priori defined constants, than the formula for
% lambda can be simplified to lamda = C * 1/tau where C = log(df/d0)
% C therefore acts as a type of physical scaling parameter (stretching)
% Requires:
% 'distMaxGroup' from particleAdvectLCS
% 'x0n' 'y0n' initial particle seed position matrix

%% INITIAL PARAMETERS
d0 = delta0;
df = 2*delta0;
C = log(df/d0);
%%

%% STEP 1: CHECK TO SEE THAT df IS NOT SET TOO LARGE FOR INPUT DATA SET
for itxdf = 1:n;
numdf(1,itxdf) = size(find(distMaxGroup(:,itxdf)>=df),1); % calculate the number of neighboring particles larger than the threshold distance per time step
fractiondf = size(numdf,1)/num0; 
end
disp(['Fraction of points >= df: ' num2str(fractiondf)]);
%%