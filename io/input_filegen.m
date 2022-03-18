clc;
clear all;
close all
format long g;

% Mesh file generation:
poly_type = 'ORTH_LEGENDRE';
% poly_type = 'BERNSTEIN';
% poly_type = 'ORTH_CHEBYSHEV';     % -- didn't work for scattered data
grid_type = 'UNSTRUCTURED';
M_arr = int64([5]);
len_M_arr = length(M_arr);
MAX_ITER = 10;

if strcmp(grid_type, 'UNSTRUCTURED') == true,
    % Generate the folllowing file before hand from available scattered
    % point data (Yudi to supply) This can be of any dimension
    XFILE = 'xR6.in'; %'g_x1.dat'
    prompt1 = 'Input the data (comma-separated, *.dat) file containing the coordinates: ';
    fileName = input(prompt1);
    coord_data = csvread(fileName);
    % NOTE: numnod = sample_size will be read from XFILE in this case. 
    [numnod, nsdim] = size(coord_data);
elseif strcmp(grid_type, 'structured') == true,
    if (nsdim == int64(1))
        switch poly_type
            case 'ORTH_CHEBYSHEV';
                %                 Chebyshev-Gauss (CG) grid.
                xx = ChebyNodes(low_lim, up_lim, double(ndivl));
                data = xx;
                numnod = length(data);
            case 'ORTH_LEGENDRE'
%                 
%               ************** WARNING ****************
%                          Uniform grid sucks
%               ************** WARNING ****************
% 
%                 xx = linspace(low_lim, up_lim, ndivl);
                xx = ChebyNodes(low_lim, up_lim, double(ndivl));
                data = xx;
                numnod = length(data);
        end
    elseif (nsdim == int64(2))
        switch poly_type
            case 'ORTH_CHEBYSHEV'
                xx = ChebyNodes(low_lim, up_lim, double(ndivl));
                yy = ChebyNodes(low_lim, up_lim, double(ndivw));
                data = cartprod(xx, yy);
                numnod = size(data,1);
            case 'ORTH_LEGENDRE'
%                 
%               ************** WARNING ****************
%                          Uniform grid sucks
%               ************** WARNING ****************                
% 
%                 xx = linspace(low_lim, up_lim, ndivl);
%                 yy = linspace(low_lim, up_lim, ndivw);                
                xx = ChebyNodes(low_lim, up_lim, double(ndivl));
                yy = ChebyNodes(low_lim, up_lim, double(ndivw));
                
                data = cartprod(xx,yy);
                numnod = ndivl*ndivw;
        end
    elseif (nsdim == int64(3))
        switch poly_type
            case 'ORTH_CHEBYSHEV'
                xx = ChebyNodes(low_lim, up_lim, double(ndivl));
                yy = ChebyNodes(low_lim, up_lim, double(ndivw));
                zz = ChebyNodes(low_lim, up_lim, double(ndivh));
                data = cartprod(xx, yy, zz);
                numnod = size(data,1);
            case 'ORTH_LEGENDRE'
% 
%               ************** WARNING ****************
%                          Uniform grid sucks
%               ************** WARNING ****************
% 
%                 xx = linspace(low_lim, up_lim, ndivl);
%                 yy = linspace(low_lim, up_lim, ndivw);
%                 zz = linspace(low_lim, up_lim, ndivh);
                xx = ChebyNodes(low_lim, up_lim, double(ndivl));
                yy = ChebyNodes(low_lim, up_lim, double(ndivw));
                zz = ChebyNodes(low_lim, up_lim, double(ndivh));
                data = cartprod(xx, yy, zz);
                numnod = ndivl*ndivw*ndivh;
        end
    end
    
end

%  Solvers for LS (all in the Expert category):
% solver = 'DGELSS';  % This one doesn't converge for Legendre and Bernstein
solver = 'DGELS';

% Solver for Normal Equation:
% solver = 'DSYSV';
% solver = 'DPOSV'; % This one does't work that well for Legendre and Bernstein
% solver = 'DGESV';
% Expert routines for Normal Equation:
% solver = 'DGESVX';

strdim = num2str(nsdim);

strn = num2str(numnod);

func_type_num = int64(10);

switch func_type_num
    case 1
        func_type_str = 'sin_x2_pl_dotdot';
    case 2
        func_type_str = 'oscillatory';
    case 3
        func_type_str = 'cos_x2_pl_dotdot';
    case 4
        x0 = 0.5*ones(1,nsdim);
        cfactors = sqrt(2.0)*ones(1,nsdim);
        alpha = 15.0;
        func_type_str = ['Gaussian_alpha_' num2str(alpha)];
    case 5
        func_type_str = 'product_peak';
    case 6
        func_type_str = 'corner_peak';
    case 7
        wave_num = 20.0;
        Q = 1;
        nn = 0;
        func_type_str = ['PU_ENRICH_COS_Q_' num2str(Q) '_DIR_' num2str(nn+1) ...
            '_Ang_' num2str(360.0*nn/Q+45) '_deg'];
    case 8
        wave_num = 20.0;
        Q = 1;
        nn = 0;
        func_type_str = ['PU_ENRICH_SIN_Q_' num2str(Q) '_DIR_' num2str(nn+1) ...
            '_Ang_' num2str(360.0*nn/Q+45) '_deg'];
        
    case 9
        l = 5;
        m = 3;
        x0 = [0.0, 0.0, 0.0];
        alpha = 0.5;
        func_type_str = ['Exp_sph_harm_alpha_' num2str(alpha) ...
                         'ELL_' num2str(l) '_EMM_' num2str(m)];
    case 10
        func_type_str = 'traffic';
%         GXFILE = 'g_x1.in';
%         GXFILE = 'g_x2.in';
%         GXFILE = 'g_x3.in';
%         GXFILE = 'g_x4.in';
        GXFILE = 'g_x5.in';
    otherwise
        disp('soon more to come!');
        error('but not just coded yet');
end
% TOL = 1.0E-03;
TOL = 1.0E-05;
C_min = -10.0;
C_max = 10.0;

% prblm_type = 'NRML_EQTN';
prblm_type = 'LLS';

if strcmpi(prblm_type, 'NRML_EQTN')
    t_reg_state = 'ON';
    lam = 100.0*eps;
elseif strcmpi(prblm_type, 'LLS')
    t_reg_state = 'OFF';
    lam = 0.0;
end

R_max = int64(6);

int_order = int64((1 + max(M_arr))/2);

COORDFILE = [strdim 'd_' strn '_Nodes_' grid_type '_mesh.dat'];

RUNDATAFILE = [strdim 'd_' strn '_Nodes_' grid_type '_mesh.runtime'];

OUTFILENAME = [strdim 'd_' strn '_Nodes_' grid_type '_mesh_' ...
    func_type_str '.out'];

PLOTFILE = [strdim 'd_' strn '_Nodes_' grid_type '_mesh_' ...
    func_type_str '_MATLABplot.dat'];

NFINEDIV = 100;

filename1 = [strdim 'd_' strn '_Nodes_' grid_type '_mesh_' ...
    func_type_str '.in'];

fileID1 = fopen(filename1,'w');
% fprintf(fileID1,'%30s \n','DIMENSION');
% fprintf(fileID1,'%d \n',nsdim);
% fprintf(fileID1,'%30s\n','SAMPLE SIZE');
% fprintf(fileID1,'%d\n',numnod);
% 
fprintf(fileID1,'%9s\n','GRID_TYPE');
fprintf(fileID1,'%30s\n', grid_type);
fprintf(fileID1,'%30s\n', 'DATA SITES');
fprintf(fileID1, '%30s\n', XFILE);
% 
fprintf(fileID1,'%30s\n','ALS MAXITER');
fprintf(fileID1,'%d\n',MAX_ITER);
% fprintf(fileID1,'%30s \n','DOMAIN');
% fprintf(fileID1,'%20.16f %20.16f\n', up_lim, low_lim);
% fprintf(fileID1,'%30s\n','COORD');
% fprintf(fileID1,[repmat('%20.16f', 1, nsdim) '\n'], data');
fprintf(fileID1,'%30s\n','FUNCTION');
fprintf(fileID1,'%d \n',func_type_num);
if func_type_num == 4,
    fprintf(fileID1,'%30s\n','CENTER OF GAUSSIAN');
    fprintf(fileID1,'%20.16f \n', x0);
    fprintf(fileID1,'%30s\n','VECTOR DECAY FACT');
    fprintf(fileID1,'%20.16f \n', cfactors);
    fprintf(fileID1,'%30s\n','SCALAR DECAY FACTOR');
    fprintf(fileID1,'%20.16f \n',alpha);
end
if (func_type_num == 7) || (func_type_num == 8)
    fprintf(fileID1,'%30s\n','WAVE NUMBER');
    fprintf(fileID1,'%20.16f \n',wave_num);
    fprintf(fileID1,'%30s\n','# OF PLANE-WAVE ENRICH.');
    fprintf(fileID1,'%d \n',Q);
    fprintf(fileID1,'%30s\n','DIR_NUM');
    fprintf(fileID1,'%d \n',nn);
end
%           
if func_type_num == 9,
    fprintf(fileID1,'%30s\n','DEG. of LEG_POLY (l)');
    fprintf(fileID1,'%d \n',l);
    fprintf(fileID1,'%30s\n','ORD. of LEG_POLY (m)');
    fprintf(fileID1,'%d \n',m);
    fprintf(fileID1,'%30s\n','CENTER OF GAUSSIAN');
    fprintf(fileID1,'%20.16f \n', x0);
    fprintf(fileID1,'%30s\n','DECAY FACTOR');
    fprintf(fileID1,'%20.16f \n',alpha);
end
if func_type_num == 10,
    fprintf(fileID1,'%30s\n','FUNC FILE NAME');
    fprintf(fileID1,'%30s\n',GXFILE);
end
% 
fprintf(fileID1,'%30s \n','ERROR TOL');
fprintf(fileID1,'%20.16f \n',TOL);
fprintf(fileID1,'%30s\n','LEN OF POLY ORDER ARRAY');
fprintf(fileID1,'%d \n',len_M_arr);
fprintf(fileID1,'%30s\n','POLY ORDER ARRAY');
fprintf(fileID1,'%d \n', M_arr);
fprintf(fileID1,'%30s \n','MAX RANK');
fprintf(fileID1,'%d \n', R_max);
fprintf(fileID1,'%30s \n','POLY TYPE');
fprintf(fileID1,'%30s \n', poly_type);
fprintf(fileID1,'%30s \n','INTEGRATION ORDER');
fprintf(fileID1,'%d \n',int_order);
fprintf(fileID1,'%30s\n','C_MIN');
fprintf(fileID1,'%20.16f \n',C_min);
fprintf(fileID1,'%30s \n','C_MAX');
fprintf(fileID1,'%20.16f \n',C_max);
fprintf(fileID1,'%30s \n','LAPACK SOLVER');
fprintf(fileID1,'%30s \n', solver);
fprintf(fileID1,'%30s \n','PROBLEM TYPE');
fprintf(fileID1,'%30s \n', prblm_type);
fprintf(fileID1,'%30s \n','TIKHONOV REG. STATE');
fprintf(fileID1,'%30s \n', t_reg_state);
fprintf(fileID1,'%30s \n','REGULARIZATION PARAM.');
fprintf(fileID1,'%20.16f \n', lam);
fprintf(fileID1,'%30s \n','OUT FILE NAME');
fprintf(fileID1,'%30s \n', OUTFILENAME);
fprintf(fileID1,'%30s \n','RUNTIME DATA FILE');
fprintf(fileID1,'%30s \n', RUNDATAFILE);
fprintf(fileID1,'%30s \n','MESH FILE');
fprintf(fileID1,'%30s \n', COORDFILE);
fprintf(fileID1,'%30s \n','PLOT FILE');
fprintf(fileID1,'%30s \n', PLOTFILE);
fprintf(fileID1,'%30s \n','RESOLUTION');
fprintf(fileID1,'%d \n',NFINEDIV);
fclose(fileID1);
