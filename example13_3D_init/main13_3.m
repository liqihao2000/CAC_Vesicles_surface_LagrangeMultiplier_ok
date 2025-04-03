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

N  = 128;
Nx = N;
Ny = N;
Nz = N;

% Parameters
para.C0 = 10000; % SAV
para.epsilon = 6*pi/128;
% para.epsilon = 0.12;
para.M = 1;
para.M1 = 50;
para.M2 = 1000;
% para.alpha = 1;
% para.alpha = 0.1;
para.alpha = 0.0001;

para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

% Time: dt T
T = 4;
t0 = 0;
tsave = 0.04;
% tsave = 0.2;

% dt_array = 0.01./2.^(1);
dt_array = 0.01./2.^(2:9);
% dt_array = 0.01./2.^(0:1);
dt_ref = 5e-4;
dt_ref = 1e-3;
dt_ref = 1e-2;

maxIt = length(dt_array);

%% option
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 1;
option.savefinal  = 0;
option.energyflag = 1;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

%% user parameters
data_array = {'ex13_Vesicles_data_eight','ex13_Vesicles_data_six',...
    'ex13_Vesicles_data_five','ex13_Vesicles_data_four',...
    'ex13_Vesicles_data_three_elliptic','ex13_Vesicles_data_threeOoo',...
    'ex13_Vesicles_data_twoOO',...
    };

% data_array = {'ex13_Vesicles_data_eight'};
index_fig = 1;
for index = 1:length(data_array)
    pdename = str2func(data_array{index});

    scheme0 = 'CAC_Vesicle_3D_';

%     scheme1_array = {'MSAV_'};
    scheme1_array = {'MSAV_','LM0_SAV_','LM1_SAV_','LM3_LM_','LM4_EX_LM_'};
%     scheme1_array = {'LM0_SAV_'};
    scheme1_array = {'LM4_EX_LM_'};

    scheme1_array = {'LM0_SAV_'};
    

    % scheme2_array = {'1st','BDF'};
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
            if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
                time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
                solver_fun(pde,domain,Nx,Ny,Nz,time,option);
            end
        end
    end
end

