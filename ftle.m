%ftle - calculates finite TIME lyapunov exponent 
disp('FTLE - calculate the finite time Lyapunov exponent.');

%% FTLE USING SHADDEN'S GRADIENT METHOD
% [ A B ; C D]
% NOTE: for most cases, the displacement from seed to secondary particle is
% fixed for up,dn,lt, and rt secondary particles, in such cases, the
% partial differentials (Jacobian) in the denominator can simply be set
% to the initial secondary perturbation distance "2*delta" in ADVECT
% The col values corespond to the END of each TS
%
%             UP
%           delta
% LT delta (xnS,ynS) delta RT
%           delta
%             DN
%
%% CALCULATE JACOBIAN USING CENTRAL DIFFERENCES
A = (xnRt(:,1:end) - xnLt(:,1:end))./2*delta; % x difference rt - lt
B = (xnUp(:,1:end) - xnDn(:,1:end))./2*delta; % x difference up - dn 
C = (ynRt(:,1:end) - ynLt(:,1:end))./2*delta; % y difference rt - lt
D = (ynUp(:,1:end) - ynDn(:,1:end))./2*delta; % y difference up - dn
%%

% NOTE: because we take the ln of the sqrt(max(eig(phiTphi)) in the FTLE
% calculation, this term must be significantly larger than e
T = abs(dt*TS); % take absolute value since FTLE can be run backwards!

for itFTLECol = 1:size(A,2);
    for itFTLE = 1:size(A,1);
        % create 2 x 2 matrix of the Jacobian from central differences
        phi = [A(itFTLE,itFTLECol) B(itFTLE,itFTLECol) ; C(itFTLE,itFTLECol) D(itFTLE,itFTLECol)];
        % set NaN values to zero
        phi(isnan(phi)) = 0;
        % left multiply the transpose of the Jacobian
        phiTphi = phi'*phi;
        % find max eigenvalue and take square root
        d(itFTLE,itFTLECol) = sqrt(max(eig(phiTphi)));
        % calculate the FTLE
        ftleout(itFTLE,itFTLECol) = log(max(eig(phiTphi)))/(T);
    end
    disp(['TS: ', num2str(itFTLECol), ' FTLE calculated.']);
end

% set Inf or -Inf values to zero
ftleout(isinf(ftleout(:))) = 0;
%%

% scale FTLE field from 0 to 1
ftleMin = min(ftleout(:));
%ftleMin = 1.0;
ftleMax = max(ftleout(:));
ftleRange = ftleMax - ftleMin;
ftleScale = 1/ftleRange;
ftleSc = (ftleout - ftleMin).*ftleScale; % transforms to a zero base
%