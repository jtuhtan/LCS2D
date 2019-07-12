% gyrefield - creates test double gyre velocity field which can be used to investigate LCS

%%
% Jeff Tuhtan 2013 - Distributable GPL
% jtuhtan@gmail.com
% V.20130208
%%

clear all
close all

filePath = 'C:\Users\tuhtan\Downloads\MATLAB\LCS\';
fileName = 'gyrefieldOUT.mat';
fileFull = strcat(filePath,fileName);

% assign x bounds
xMin = 0;
xMax = 2;
% assign y bounds
yMin = 0;
yMax = 1;

resGridVel = 0.005; % velocity field resolution, same in x and y directions

[gX, gY] = meshgrid(xMin:resGridVel:xMax,yMin:resGridVel:yMax); % create grid 

psi = sin(pi*gX).*sin(pi*gY); % create psi - scalar potential field

xyz = [gX(:),gY(:),psi(:)]; % create column vector format of points [x,y,psi]

[gradX, gradY] = gradient(psi); % calculate gradient of scalar potential field

u = -gradY; % retrieve u velocity component (double gyre eqn.)
v = gradX; % retrive v velocity component (double gyre eqn.)

velMag = bsxfun(@hypot,u,v); % calculate velocity magnitude 2D

velocity = [u(:),v(:)]; % assign to column array for later use

xyzvelocity = [xyz,velocity]; % assign to xyzuv column array 

%% Plot potential surface and velocity vectors
figure(1);
contourf(gX,gY,velMag), hold on, quiver(gX,gY,u,v), axis equal, axis tight, 
xlabel('x [m]');
ylabel('y [m]');
text(2.48,0.9,'Velocity Magnitude [m/s]','rot',-90);
caxis([0, 0.08]);
cbar = colorbar;

figure(2);
contourf(gX,gY,u), hold on, quiver(gX,gY,u,v), axis equal, axis tight, 
xlabel('x [m]');
ylabel('y [m]');
text(2.50,0.78,'X Velocity [m/s]','rot',-90);
caxis([-0.1, 0.1]);
cbar = colorbar;

figure(3);
contourf(gX,gY,v), hold on, quiver(gX,gY,u,v), axis equal, axis tight, 
xlabel('x [m]');
ylabel('y [m]');
text(2.50,0.78,'Y Velocity [m/s]','rot',-90);
caxis([-0.1, 0.1]);
cbar = colorbar;
%%

save fileFull

%% OPTIONAL: Export data file containing gyre nodal locations
% At home
%dlmwrite('K:\LCS\gyrenodes.xyz', xyz, 'delimiter', '\t', 'precision', 6);
%dlmwrite('K:\LCS\gyrevelocity.xyz', velocity, 'delimiter', '\t', 'precision', 6);         
%%