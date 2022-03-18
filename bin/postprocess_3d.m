%  Some plotting:
clc;
clear all;
% close all;
format long g;

name = '3d_421875_Nodes_structured_mesh_Exp_sph_harm_alpha_0.5ELL_5_EMM_3_MATLABplot.dat';
fid = fopen(name, 'r');
formatSpec = '%f %f %f %f';
sizeofarray = [4 inf^2];
matrixdata = fscanf(fid, formatSpec, sizeofarray);
matrixdata = matrixdata';

XX = matrixdata(:,1);
YY = matrixdata(:,2);
ZZ = matrixdata(:,3);
VAL = matrixdata(:,4);

xvec = unique(XX);
yvec = unique(YY);
zvec = unique(ZZ);

npts = int64((length(XX))^(1/3));

% X = reshape(XX, [npts, npts, npts]);
% Y = reshape(YY, [npts, npts, npts]);
% Z = reshape(ZZ, [npts, npts, npts]);

[X, Y, Z] = meshgrid(xvec, yvec, zvec);

R_times_Y = reshape(VAL, [npts, npts, npts]);

alw = 0.75;    % AxesLineWidth
fsz = 25;      % Fontsize
lw = 1.5;      % LineWidth

RY_max = max(R_times_Y(:));     RY_min = min(R_times_Y(:));

figure
clf;
vals_I = linspace(RY_max, RY_min, 50);     % 20 points from max down to min
h_I = patch(isosurface(X, Y, Z, R_times_Y, vals_I(1)), ...
              'facecolor', 'g', 'edgecolor', 'none');

axis equal
xlim([min(min(min(X))) max(max(max(X)))]);
ylim([min(min(min(Y))) max(max(max(Y)))]);
zlim([min(min(min(Z))) max(max(max(Z)))]);
xlim([min(min(min(X))) max(max(max(X)))]);
ylim([min(min(min(Y))) max(max(max(Y)))]);
zlim([min(min(min(Z))) max(max(max(Z)))]);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'zlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');
axis vis3d
axis normal
axis equal
light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[10 -1 10]);
view (40,34);
grid on
lighting phong
camlight('headlight');

for id = vals_I
   [faces, vertices] = isosurface(X, Y, Z, R_times_Y, id);
   set(h_I, 'Faces', faces, 'Vertices', vertices);
   pause(0.1);
end



% % Plotting
% width = 4;     % Width in inches
% height = 2;    % Height in inches
% alw = 0.75;    % AxesLineWidth
% fsz = 25;      % Fontsize
% lw = 1.5;      % LineWidth
% msz = 8;       % MarkerSize
% 
% figure;
% set(gcf,'InvertHardcopy','on');
% set(gcf,'PaperUnits', 'inches');
% pos = get(gcf, 'Position');
% set(gca, 'FontSize', fsz, 'LineWidth', alw);
% 
% isosurface(X,Y,Z,fvals);
% 
% xlim([min(min(min(X))) max(max(max(X)))]);
% ylim([min(min(min(Y))) max(max(max(Y)))]);
% zlim([min(min(min(Z))) max(max(max(Z)))]);
% set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% set(get(gca,'zlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% set(gca,'LineWidth',lw);
% set(gca,'FontSize',fsz);
% set(gca,'FontWeight','Bold');
% set(gcf,'color','w');
% axis vis3d
% axis normal
% axis equal
% light('Position',[-50 -15 29]);
% set(gca,'CameraPosition',[10 -1 10]);
% view (40,34);
% lighting phong
% grid on
% camlight('headlight');
% % filename_enr = ['Rn' num2str(l) '_times_Y' num2str(m) num2str(l)];
% % print(filename_enr,'-depsc2');
% 
% 
% 
% % % Little triangles
% % % The solution is to use Delaunay triangulation. Let's look at some
% % % info about the "tri" variable.
% % tri = delaunay(XX',YY');
% % figure
% % plot(XX',YY','r.','MarkerSize', 1)
% % hold on
% % % % How many triangles are there?
% % % [r,c] = size(tri);
% % % disp('# of traiangles resulting from Delaunay: ')
% % % disp(r)
% % 
% % % Plot it with TRISURF
% % h = trisurf(tri, XX', YY', ZZ');
% % axis vis3d
% % % Clean it up:
% % axis normal
% % % axis off
% % l = light('Position',[-50 -15 29]);
% % set(gca,'CameraPosition',[10 -1 10]);
% % view (50,30);
% % lighting phong
% % shading interp
% % colorbar EastOutside
% % camlight
% % grid on
% % light;
% % title('Separated Representation (Surface Plot)','FontSize',14);
% % xlabel('X')
% % ylabel('Y')
% % zlabel('Z')
% % % goodplot
% % 
% % % wvn = 20.0;
% % % qq = 1;
% % % nn = 0;
% % % theta_n = 2.0*pi*nn/qq + pi/4;
% % % % f7 = @(x,y) cos(wvn*(x.*cos(theta_n) + y.*sin(theta_n)));
% % % % f_ana = f7(XX,YY);
% % % f8 = @(x,y) sin(wvn*(x.*cos(theta_n) + y.*sin(theta_n)));
% % % f_ana = f8(XX,YY);
% % 
% % t = 0.5;
% % C = sqrt(2.)*ones(2,1);
% % w = [1./8.; 1./8.];
% % f4 = @(x,y) exp(-t*(C(1)^2*(x - w(1)).^2 + C(2)^2*(y - w(2)).^2));
% % f_ana = f4(XX,YY);
% % 
% % 
% % figure
% % plot(XX',YY','r.','MarkerSize', 1)
% % hold on
% % 
% % hh = trisurf(tri, XX', YY', f_ana');
% % axis vis3d
% % % Clean it up:
% % axis normal
% % % axis off
% % l = light('Position',[-50 -15 29]);
% % set(gca,'CameraPosition',[10 -1 10]);
% % view (50,30);
% % lighting phong
% % shading interp
% % colorbar EastOutside
% % camlight
% % grid on
% % light;
% % title('Analytical Representation of the function (Surface Plot)','FontSize',14)
% % % goodplot
% % % res_Helm2d = integral2(f7, -1.0, 1.0, -1.0, 1.0);
% % % num = 1.9764401586253477E-002;
% % res_Helm2d = integral2(f4, -1.0, 1.0, -1.0, 1.0);
% % num = 2.1984070448510806E-005;
% % err_num = abs(res_Helm2d - num)/abs(res_Helm2d)*100.0;
% % fprintf('The Error in integration: %20.16f%% \n', err_num)