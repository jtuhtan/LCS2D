% ftleplot - plots view of FTLE 
disp('FTLEPLOT - visualize the birth of the FTLE before your own eyes.');
close all

%% CONTOUR PLOT OF FTLE FIELD
% open video and set frame rate and quality
aviobj = VideoWriter('C:\Users\tuhtan\Downloads\MATLAB\LCS\Videos\ftleplot.mp4');
aviobj.FrameRate = 30;
aviobj.Quality = 100;
% open video file for writing
open(aviobj);

ptStart = 3*size(seedX,1)+1;
ptStop = ptStart+size(seedX,1)-1;
plotInt = 1; % plot interval for high-density seeds
MarkerSz = 2; % markersize in plot

idxSeedNaN = find(~isnan(xnS(:,1))); % index filtered seeds
f1 = figure(1);
set(f1,'Position',[50 450 750 450]);

for itFTLECon = 2:TS;
ftleF = TriScatteredInterp(xnS(idxSeedNaN,1),ynS(idxSeedNaN,1),ftleout(idxSeedNaN,itFTLECon),'natural'); % velocity and x y datasets must correspond to the same nodes!
ftleG(:,:,itFTLECon) = ftleF(seedX,seedY);

now = ftleG(:,:,itFTLECon);
if itFTLECon > 1;
previous = ftleG(:,:,itFTLECon-1);
idx = find(now==0);
now(idx) = previous(idx);
else
end
ftleG(:,:,itFTLECon) = now;

contourf(seedX,seedY,now,30,'LineStyle','none');
%surf(seedX,seedY,now,'FaceColor','interp','EdgeColor','None','FaceLighting','Phong');
%surfl(seedX,seedY,now);
%shading interp
%hold on;

%contourf(seedX,seedY,magSeed,10,'edgecolor','none');
%surf(seedX,seedY,magSeed./10000);
colormap jet
%shading interp
hold on

%h1 = plot3(xnS(1:plotInt:end,itFTLECon),ynS(1:plotInt:end,itFTLECon),0*ynS(1:plotInt:end,itFTLECon)+1,'o');
%set(h1, 'Markersize',MarkerSz);
%set(h1,'MarkerFaceColor',[0.0 0.0 0.0]);
%set(h1,'MarkerEdgeColor','none');
%hold on % turn on to get 'tracer effect'
plotSecondary = 0;
if plotSecondary == 1;
    h2 = plot3(xnS(1:plotInt:end,itFTLECon),ynS(1:plotInt:end,itFTLECon),D(1:plotInt:end,itFTLECon),'.c');
    set(h2, 'Markersize',MarkerSz);
    hold on
    h3 = plot3(xnUp(1:plotInt:end,itFTLECon),ynUp(1:plotInt:end,itFTLECon),0*ynS(1:plotInt:end,itFTLECon)+2,'.r');%
    set(h3, 'Markersize',MarkerSz);
    h4 = plot3(xnDn(1:plotInt:end,itFTLECon),ynDn(1:plotInt:end,itFTLECon),0*ynS(1:plotInt:end,itFTLECon)+2,'.g');
    set(h4, 'Markersize',MarkerSz);
    h5 = plot3(xnLt(1:plotInt:end,itFTLECon),ynLt(1:plotInt:end,itFTLECon),0*ynS(1:plotInt:end,itFTLECon)+3,'.y');
    set(h5, 'Markersize',MarkerSz);
    h6 = plot3(xnRt(1:plotInt:end,itFTLECon),ynRt(1:plotInt:end,itFTLECon),0*ynS(1:plotInt:end,itFTLECon)+4,'.b');
    set(h6, 'Markersize',MarkerSz);
else
end
%colorbar
axis equal
set(gca,'Color','k');
set(gcf,'Color',[0 0 0]);
%camlight right
view(2);
%axis equal
%caxis([-70, -60]);
%zlim([-35 -34]);
%xlim([0,1]);
%ylim([0,1]);
pause(0.001);
hold off

% save movie
frame = getframe(gcf);
writeVideo(aviobj,frame);

clear ftleF
end
close(aviobj);
%%

