close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
domain.xleft   = 0;
domain.xright  = 2*pi;
domain.yleft   = 0;
domain.yright  = 2*pi;
domain.zleft   = 0;
domain.zright  = 2*pi;

Lx = domain.xright - domain.xleft;
Ly = domain.yright - domain.yleft;
Lz = domain.zright - domain.zleft;

N  = 64;
Nx = N;
Ny = N;
Nz = N;

% Parameters
para.C0 = 1000; % SAV
para.epsilon = 6*pi/128;
% para.epsilon = 0.07;
para.M = 1;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

% Time: dt T
T = 2;
t0 = 0;
tsave = 0.5*T;

% dt_array = 0.01./2.^(1); 
dt_array = 0.01./2.^(2:6); 
% dt_array = 0.01./2.^(0:1);
dt_ref = 4e-4;

maxIt = length(dt_array);

%% option
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 0;
option.savefinal  = 0;
option.energyflag = 1;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

%% user parameters
pdename = str2func('ex11_Vesicles_data');
scheme0 = 'CAC_Vesicle_3D_';

scheme1_array = {'LM0_SAV_','LM1_SAV_','LM3_LM_'};
scheme1_array = {'LM3_LM_'};

scheme2_array = {'1st','BDF'};
scheme2_array = {'BDF'};

for i_2 = 1:length(scheme2_array)
    for i_1 = 1:length(scheme1_array)
    
        scheme1 = scheme1_array{i_1};
        scheme2 = scheme2_array{i_2};

        scheme = [scheme0,scheme1,scheme2]
        para.name = [func2str(pdename),'_',[scheme1,scheme2]];

        pde = pdename(para);
        solver_fun = str2func(scheme);
        %% Run:
        for k = 1:maxIt
            if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
                time = struct('T',T,'t0',t0,'dt',dt_array(k),'tsave',tsave);
                solver_fun(pde,domain,Nx,Ny,Nz,time,option);
            end
        end
    end
end

