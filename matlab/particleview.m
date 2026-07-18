% particleview - view PAM results
disp('PARTICLEVIEW - view advected particles.');

%close all
pointSize = 5;

for itTS = 1:TS;
%imagesc(A);
%% Solid bg color
whitebg('k');
%%
colormap gray
%hold on % activate this to see the integrated motion from each time step
h1 = plot(xnS(:,itTS),ynS(:,itTS),'.c');
set(h1,'MarkerSize',pointSize);
hold on
h2 = plot(xnUp(:,itTS),ynUp(:,itTS),'.r');
set(h2,'MarkerSize',pointSize);
h3 = plot(xnDn(:,itTS),ynDn(:,itTS),'.g');
set(h3,'MarkerSize',pointSize);
h4 = plot(xnLt(:,itTS),ynLt(:,itTS),'.b');
set(h4,'MarkerSize',pointSize);
h5 = plot(xnRt(:,itTS),ynRt(:,itTS),'.c');
set(h5,'MarkerSize',pointSize);
axis([min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2))]);
axis equal
pause(0.01);
hold off
end

%% PLOT 
%plot(xn.A(:,TS),yn.A(:,TS),'.b');
%axis equal
%hold on
%plot(xn.B(:,TS),yn.B(:,TS),'.r');
%%