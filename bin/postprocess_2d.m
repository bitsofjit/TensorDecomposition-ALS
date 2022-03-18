%  Some plotting:
clc;
clear all;
close all;
format long g;
name = '2d_10000_Nodes_structured_mesh_Gaussian_alpha_15_MATLABplot.dat';
fid = fopen(name, 'r');
formatSpec = '%f %f %f';
sizeofarray = [3 inf];
matrixdata = fscanf(fid, formatSpec, sizeofarray);
matrixdata = matrixdata';

XX = matrixdata(:,1);
YY = matrixdata(:,2);
ZZ = matrixdata(:,3);

% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.
tri = delaunay(XX',YY');
figure
plot(XX',YY','r.','MarkerSize', 1)
hold on
% % How many triangles are there?
% [r,c] = size(tri);
% disp('# of traiangles resulting from Delaunay: ')
% disp(r)

% Plot it with TRISURF
h = trisurf(tri, XX', YY', ZZ');
axis vis3d
% Clean it up:
axis normal
% axis off
l = light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[10 -1 10]);
view (50,30);
lighting phong
shading interp
colorbar EastOutside
camlight
grid on
light;
title('Separated Representation (Surface Plot)','FontSize',14);
xlabel('X')
ylabel('Y')
zlabel('Z')
% goodplot

% wvn = 20.0;
% qq = 1;
% nn = 0;
% theta_n = 2.0*pi*nn/qq + pi/4;
% % f7 = @(x,y) cos(wvn*(x.*cos(theta_n) + y.*sin(theta_n)));
% % f_ana = f7(XX,YY);
% f8 = @(x,y) sin(wvn*(x.*cos(theta_n) + y.*sin(theta_n)));
% f_ana = f8(XX,YY);
w = 0.5*ones(2,1);
C = sqrt(2.0)*ones(2,1);
t = 15.0;
f4 = @(x,y) exp(-t*(C(1)^2*(x - w(1)).^2 + C(2)^2*(y - w(2)).^2));
f_ana = f4(XX,YY);


figure
plot(XX',YY','r.','MarkerSize', 1)
hold on

hh = trisurf(tri, XX', YY', f_ana');
axis vis3d
% Clean it up:
axis normal
% axis off
l = light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[10 -1 10]);
view (50,30);
lighting phong
shading interp
colorbar EastOutside
camlight
grid on
light;
title('Analytical Representation of the function (Surface Plot)','FontSize',14)
% goodplot
% res_Helm2d = integral2(f7, -1.0, 1.0, -1.0, 1.0);
% num = 1.9764401586253477E-002;
% res_Helm2d = integral2(f4, -1.0, 1.0, -1.0, 1.0);
% num = 2.1984070448510806E-005;
% err_num = abs(res_Helm2d - num)/abs(res_Helm2d)*100.0;
% fprintf('The Error in integration: %20.16f%% \n', err_num)