% gertplot - view Gert's data quick!

%% CONTOUR PLOT OF FTLE FIELD
makeVideo = 0;
if makeVideo == 1;
% open video and set frame rate and quality
aviobj = VideoWriter('C:\Users\tuhtan\Downloads\MATLAB\LCS\Videos\gertplotVelSlow.mp4');
aviobj.FrameRate = 10;
aviobj.Quality = 100;
% open video file for writing
open(aviobj);
else
end

% Load Gert Toming PIV data, no z value!
fileName = 'C:\Users\tuhtan\Downloads\MATLAB\LCS\Data\PIV\behind the fish tail 30cms_10cm_cylinder\data\vel_field';
%fileName = 'D:\2012_12_03_Email_G_Toming_Server\KVS_10cyl_30cms_20dummy\vel_field';
%fileName = 'D:\2012_12_03_Email_G_Toming_Server\UF_30cms_20dummy\vel_field';
%fileName = 'D:\2012_12_03_Email_G_Toming_Server\UF_30cms_48dummy\vel_field';
for it = 1:1;
    load(strcat(fileName,num2str(it),'.mat'));
    disp(['vel_field',num2str(it)]);
%load D:\2012_12_03_Email_G_Toming_Server\KVS_10cyl_30cms_20dummy\vel_field3.mat
velocity = [fu(:)./100, fv(:)./100];
%velocity(isnan(velocity)) = 0;

[vort, angvel] = curl(fu./100,fv./100);
%vort(isnan(vort)) = 0;

numX = size(u,1);
numY = size(u,2);
[xGrid, yGrid] = meshgrid(1:79,1:63);
xGrid = flipud(xGrid);
yGrid = flipud(yGrid);
xyz = [xGrid(:), yGrid(:), 0*yGrid(:)];
vMag = sqrt(fu.^2 + fv.^2)./100;

if makeVideo == 1;
h = figure();
else
%figure();
surf(xGrid,yGrid,vMag,'EdgeColor','none');
%shading interp
set(gca,'Color','k');
set(gcf,'Color',[0 0 0]);
colorbar
axis equal
%caxis([-0.2 0.2]);
caxis([0 0.5]);
view(2);
pause(0.1);
end
if makeVideo == 1;
writeVideo(aviobj,getframe(h));
else
end
end
if makeVideo == 1;
close(aviobj);
else
end