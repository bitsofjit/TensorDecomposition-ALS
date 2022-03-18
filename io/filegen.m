clc; clear all; close all;

prompt1 = 'Input the data file name: ';
fileName = input(prompt1);
g = csvread(fileName);
% 
dim = size(g,2);
% 
for i = 1:dim
    % Alternatively:
    % fileID = fopen(['g_x' num2str(i) '.dat'], 'w');
    fileID = fopen(sprintf('g_x%i.in', i), 'w');
    fprintf(fileID,'%4.4f\n', g(:,i));
    fclose(fileID);
end



% prompt1 = 'Input the data file name: ';
% fileName = input(prompt1);
% x = csvread(fileName);
% % 
% nSample = size(x,1);
% % 
% fileID = fopen('xR6.in', 'w');
% 
% for i = 1:nSample
%     % Alternatively:
%     % fileID = fopen(['g_x' num2str(i) '.dat'], 'w');
%     % fileID = fopen(sprintf('g_x%i.dat', i), 'w');
%     fprintf(fileID,'%5.5f %5.5f %5.5f %5.5f %5.5f %5.5f\n', x(i,:));
%     
% end
%     
% fclose(fileID);
    
    


    
    
    

