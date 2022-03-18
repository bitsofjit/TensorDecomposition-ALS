% function plotEnrichSchEq(N_1d, xup, yup, zup, xlow, ylow, zlow, ...
%                          l, m, x0, alpha)
% 
clc;
% clear all;
% if nargin < 11,
    N_1d = 50;
%     xup = 10.0;  yup = 10.0;  zup = 10.0;
%     xlow = -10.0;  ylow = -10.0;  zlow = -10.0;                     
    xup = 1.0;  yup = 1.0;  zup = 1.0;
    xlow = 0.0;  ylow = 0.0;  zlow = 0.0; 
    l = 5;
    m = 3;
    x0 = [0.0; 0.0; 0.0];
    alpha = 0.5;
% end

x = linspace(xlow, xup, N_1d);
y = linspace(ylow, yup, N_1d);
z = linspace(zlow, zup, N_1d);

[X, Y, Z] = meshgrid(x, y, z);


% x = 1:3;
% y = 1:4;
% [X,Y] = meshgrid(x,y)
% =>

% X =  [1, 2, 3;
%       1, 2, 3;
%       1, 2, 3;
%       1, 2, 3];
%   
% Y =  [1, 1, 1;
%       2, 2, 2;
%       3, 3, 3;
%       4, 4, 4];

% x = 1:3;
% y = 1:4;
% z = 1:5;
% [XX, YY, ZZ] = meshgrid(x,y,z);
% XX => [4 x 3 x 5]
% YY => [4 x 3 x 5]
% ZZ => [4 x 3 x 5]
% 

%     + + +
%    + + + | 
%   + + +  |
%  + + +   / 
% + + +   /
% + + +  /   
% + + + /
% + + +/

R_times_Y = zeros(N_1d, N_1d, N_1d);
Yml = zeros(N_1d, N_1d, N_1d);

for i = 1:N_1d,
    for j = 1:N_1d,
        for k = 1:N_1d,
            xsample = [x(j); y(i); z(k)];
            [R_times_Y(j, i, k), Yml(j, i, k)] = ...
             Rnl_times_Yml(l, m, x0, xsample, alpha);
            
        end
    end
end

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
% title_string_II = ['(Y^m_l)^2 for l =  ' num2str(l) '  and  m = ' num2str(m)];
% title(title_string_II)
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

figure;
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gca, 'FontSize', fsz, 'LineWidth', alw);
clf;
h_II= slice(X, Y, Z, R_times_Y, [], [], vals_I(1));
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
colorbar;
for id = vals_I
    delete(h_II);
    h_II= slice(X, Y, Z, R_times_Y, [], [], id);
    xlim([min(min(min(X))) max(max(max(X)))]);
    ylim([min(min(min(Y))) max(max(max(Y)))]);
    zlim([min(min(min(Z))) max(max(max(Z)))]);
    xlim([min(min(min(X))) max(max(max(X)))]);
    ylim([min(min(min(Y))) max(max(max(Y)))]);
    zlim([min(min(min(Z))) max(max(max(Z)))]);
    colorbar;
    pause(0.1);
end

% 
% 
% 
% % figure;
% % set(gcf,'InvertHardcopy','on');
% % set(gcf,'PaperUnits', 'inches');
% % set(gca, 'FontSize', fsz, 'LineWidth', alw);
% % 
% % isosurface(X,Y,Z,R_times_Y);
% % 
% % xlim([min(min(min(X))) max(max(max(X)))]);
% % ylim([min(min(min(Y))) max(max(max(Y)))]);
% % zlim([min(min(min(Z))) max(max(max(Z)))]);
% % set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% % set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% % set(get(gca,'zlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% % title_string_temp1 = ['e^{-\alpha (r - r_0)^2} Y^m_l for l =  ' num2str(l) ', m =  ' num2str(m)];
% % title_string_temp2 = [', \alpha = ' num2str(alpha) ', and r_0 = || ' num2str(x0') '||_{2}' ];
% % title_string_I = [title_string_temp1 title_string_temp2];
% % title(title_string_I)
% % set(gca,'LineWidth',lw);
% % set(gca,'FontSize',fsz);
% % set(gca,'FontWeight','Bold');
% % set(gcf,'color','w');
% % axis vis3d
% % axis normal
% % axis equal
% % light('Position',[-50 -15 29]);
% % set(gca,'CameraPosition',[10 -1 10]);
% % view (40,34);
% % lighting phong
% % grid on
% % camlight('headlight');
% % filename_enr = ['Rn' num2str(l) '_times_Y' num2str(m) num2str(l)];
% % print(filename_enr,'-depsc2');
% 
% % 
% % Use of scatter to get an idea about the values
% % figure
% % scatter3(X(:), Y(:), Z(:), 72, (Yml(:)).^2, 'filled');
% % view(-15, 35);
% % colorbar;
% 
% % figure
% % slice(X, Y, Z, Yml.^2, [], [], 0);
% % colorbar;
% 
% 
% Use of slice function and then animating:
figure;
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
set(gca, 'FontSize', fsz, 'LineWidth', alw);
clf;
vals = 0.309467908419468:-0.01:0.01;
h_III= slice(X, Y, Z, Yml.^2, [], [], vals(1));
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
colorbar;
for id = vals
    delete(h_III);
    h_III= slice(X, Y, Z, Yml.^2, [], [], id);
    xlim([min(min(min(X))) max(max(max(X)))]);
    ylim([min(min(min(Y))) max(max(max(Y)))]);
    zlim([min(min(min(Z))) max(max(max(Z)))]);
    xlim([min(min(min(X))) max(max(max(X)))]);
    ylim([min(min(min(Y))) max(max(max(Y)))]);
    zlim([min(min(min(Z))) max(max(max(Z)))]);
    colorbar;
    pause(0.1);
end
% 
% 
% 
% 
% % figure;
% % set(gcf,'InvertHardcopy','on');
% % set(gcf,'PaperUnits', 'inches');
% % set(gca, 'FontSize', fsz, 'LineWidth', alw);
% % 
% % isosurface(X,Y,Z,Yml.^2);
% % 
% % axis equal
% % xlim([min(min(min(X))) max(max(max(X)))]);
% % ylim([min(min(min(Y))) max(max(max(Y)))]);
% % zlim([min(min(min(Z))) max(max(max(Z)))]);
% % xlim([min(min(min(X))) max(max(max(X)))]);
% % ylim([min(min(min(Y))) max(max(max(Y)))]);
% % zlim([min(min(min(Z))) max(max(max(Z)))]);
% % set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% % set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% % set(get(gca,'zlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% % title_string_II = ['(Y^m_l)^2 for l =  ' num2str(l) '  and  m = ' num2str(m)];
% % title(title_string_II)
% % set(gca,'LineWidth',lw);
% % set(gca,'FontSize',fsz);
% % set(gca,'FontWeight','Bold');
% % set(gcf,'color','w');
% % axis vis3d
% % axis normal
% % axis equal
% % light('Position',[-50 -15 29]);
% % set(gca,'CameraPosition',[10 -1 10]);
% % view (40,34);
% % lighting phong
% % grid on
% % camlight('headlight');
% 
% % filename_Yml = ['Y' num2str(m) num2str(l)];
% % print(filename_Yml,'-depsc2');
% 
% 
% % maxYmlsq = max(Yml(:).^2);
% % minYmlsq = min(Yml(:).^2);
% % minYmlsq = 0.0;
% 
figure
clf;
vals = logspace(-1, -5, 25);     % 20 points from max down to min
h_IV = patch(isosurface(X, Y, Z, Yml.^2, vals(1)), ...
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
% title_string_II = ['(Y^m_l)^2 for l =  ' num2str(l) '  and  m = ' num2str(m)];
% title(title_string_II)
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


for id = vals
   [faces, vertices] = isosurface(X, Y, Z, Yml.^2, id);
   set(h_IV, 'Faces', faces, 'Vertices', vertices);
   pause(0.1);
end

% 






