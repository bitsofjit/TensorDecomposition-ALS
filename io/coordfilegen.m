clc; clear all; close all;

prompt1 = 'Input the data file containing the coordinates: ';
fileName = input(prompt1);
x = csvread(fileName);
% 
[nSample, dim] = size(x);

% 
fmtString = repmat('%5.5f ', 1, dim);
fileID = fopen(['xR' num2str(dim) '.in'], 'w');
% 
fprintf(fileID,'%11s\n','SAMPLE SIZE');
fprintf(fileID,'%d\n',nSample);
fprintf(fileID,'%9s\n','DIMENSION');
fprintf(fileID,'%d\n',dim);
fprintf(fileID,'%6s\n','COORDS');
% 
for i = 1:nSample
    % Alternatively:
    % fileID = fopen(['g_x' num2str(i) '.dat'], 'w');
    % fileID = fopen(sprintf('g_x%i.dat', i), 'w');
    fprintf(fileID, [fmtString, '\n'], x(i,:));
    %fprintf(fileID,'%5.5f %5.5f %5.5f %5.5f %5.5f %5.5f\n', x(i,:));
    
end
    
fclose(fileID);