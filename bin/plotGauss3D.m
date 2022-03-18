function plotGauss3D(N_1d, xup, yup, zup, xlow, ylow, zlow, ...
                   c, t, x0)

if nargin < 11,
    N_1d = 75;                    
    xup = 1.0;  yup = 1.0;  zup = 1.0;
    xlow = 0.0;  ylow = 0.0;  zlow = 0.0;
    t = 0.5; 
    c = sqrt(2.0)*ones(3,1);
    x0 = 1.0/8.0*ones(3,1);
end

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

exp_func = zeros(N_1d, N_1d, N_1d);

for i = 1:N_1d,
    for j = 1:N_1d,
        for k = 1:N_1d,
            xsample = [x(j); y(i); z(k)];
            exp_func(j, i, k) = exp(-2.0*t*sum(c.^2.*(xsample - x0).^2));
        end
    end
end

width = 4;     % Width in inches
height = 2;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 25;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize


figure;
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gca, 'FontSize', fsz, 'LineWidth', alw);

isosurface(X,Y,Z,exp_func);

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
lighting phong
grid on
camlight('headlight');
% filename_enr = ['Rn' num2str(l) '_times_Y' num2str(m) num2str(l)];
% print(filename_enr,'-depsc2');




