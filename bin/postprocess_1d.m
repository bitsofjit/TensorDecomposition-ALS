%  Some plotting:
clc;
clear all;
close all;
format long g;
name = '1d_100_Nodes_structured_mesh_Helmholtz_MATLABplot.dat';

fid = fopen(name, 'r');
formatSpec = '%f %f';
sizeofarray = [2 inf];
matrixdata = fscanf(fid, formatSpec, sizeofarray);

matrixdata = matrixdata';

XX = matrixdata(:,1);
ZZ = matrixdata(:,2);

% Separated Representation
plot(XX', ZZ', 'r-');
grid on
title('Separated Representation (Surface Plot)','FontSize',14);
xlabel('X')
ylabel('f(X) [Separated Representation]')
% goodplot

% Analytical Representation
wvn = 5.0;
f7 = @(x) cos(wvn.*x./sqrt(2.0));
f_ana = f7(XX);
figure
plot(XX',f_ana','k-.')
hold on
grid on;
title('Analytical Representation of the function (Surface Plot)','FontSize',14)
xlabel('X')
ylabel('f(X) [Analytical]')
% goodplot

% Numerical and Analytical Integration:
wvn = 5.0;
fun_Helm = @(x) cos(wvn*x/sqrt(2.0));
res_Helm1d = integral(fun_Helm, -1.0, 1.0);
num = -0.10856373156178745;
err_num = abs(res_Helm1d - num)/abs(res_Helm1d)*100.0



