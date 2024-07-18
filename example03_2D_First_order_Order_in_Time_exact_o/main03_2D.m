close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
% domain.left   = -pi;
% domain.right  =  pi;
% domain.bottom = -pi;
% domain.top    =  pi;

domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N  = 64;
Nx = N;
Ny = N;

% Parameters
para.Lx = Lx;
para.Ly = Ly;
para.epsilon = 6*pi/128;
para.M = 2;
para.M1 = 50;
para.M2 = 50;

delete *.mat
for kk = 1:2
    if 1 == kk
        para.S1 = 0;
        para.S2 = 0;
        para.S3 = 0;
    elseif 2 == kk
        para.S1 = 4;
        para.S2 = 4;
        para.S3 = 1;
    end

    PDE = 'data2';

    if 1 == strcmp(PDE,'data1')
        T = 0.002;
        dt_array = 0.0001./2.^(0:6)';
        dt_ref = 1e-7;
        para.C0 = 1000; % SAV
        pde = ex03_1_Vesicles_data(para);
    end

    if 1 == strcmp(PDE,'data2')
        T = 0.1;
        %     dt_array = 0.001./2.^(0:6)';
        dt_array = 0.005./2.^(0:6)';
        dt_ref = 1e-8;
        para.C0 = 1; % SAV
        pde = ex03_2_Vesicles_data(para);
    end

    % Time: dt T

    t0 = 0;
    tsave = 1*T;

    maxIt = length(dt_array);

    %% option
    option.plotflag  = 0;
    option.printflag = 0;
    option.vtkflag  = 0;
    option.saveflag  = 0;
    option.savefinal  = 1;
    option.energyflag = 0;
    option.tol = 1e-14;
    option.tolit = 1e-11;
    option.maxit = 2000;

    %% Run:
    % delete *.mat
    if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
        time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
        %         CAC_Vesicle_2D_MSAV_1st(pde,domain,Nx,Ny,time,option);
        CAC_Vesicle_2D_LM0_SAV_1st(pde,domain,Nx,Ny,time,option);
        %     CAC_Vesicle_2D_LM1_SAV_1st(pde,domain,Nx,Ny,time,option);
        %     CAC_Vesicle_2D_LM3_LM_1st(pde,domain,Nx,Ny,time,option);
    end
    for k = 1:maxIt
        dt = dt_array(k);
        time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
        %         v2 = CAC_Vesicle_2D_MSAV_1st(pde,domain,Nx,Ny,time,option);
        v2 = CAC_Vesicle_2D_LM0_SAV_1st(pde,domain,Nx,Ny,time,option);
        %     v2 = CAC_Vesicle_2D_LM1_SAV_1st(pde,domain,Nx,Ny,time,option);
        %     v2 = CAC_Vesicle_2D_LM3_LM_1st(pde,domain,Nx,Ny,time,option);
    end

    %% Compute order of convergence
    error=zeros(maxIt,1);
    order=zeros(maxIt,1);
    if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
        name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),'S1=',num2str(pde.S1),...
            'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
        filename=[name '.mat'];
        load(filename,'phi');
        phi_exact = phi;
        clear phi;
        fprintf('Accuracy test with reference solution.\n');
    end
    for k = 1:maxIt
        dt = dt_array(k);
        name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),'S1=',num2str(pde.S1),...
            'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
        filenamek=[name '.mat'];
        load(filenamek,'phi','xx','yy','hx','hy');

        if ~exist('phi_exact','var') && ~exist('psi_exact','var')
            phi_exact = pde.exact(xx,yy,T);
            fprintf('Accuracy test with exact solution.\n');
        end

        err    = fft2((phi_exact - phi).^2);
        error(k,1) = sqrt(err(1,1)*hx*hy);   % L2
        clear phi;
    end
    order(2:maxIt) = log(error(1:maxIt-1)./error(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

    %% Display error and order
    if 0 ~= para.S1
        fprintf('Stabilized method\n');
    end
    fprintf('    dt     &   Error_L2   & Order \n');
    for k = 1:maxIt
        fprintf('%.4e & %.4e & %.2f \n',dt_array(k),error(k),order(k));
        %         fprintf('1/%d\t& %.6e\t& %.4f %s \n',1./dt_array(k),error(k),order(k),'\\');
    end
    fprintf('\n')

    %% Save error and order
    name=['phi_e',num2str(para.epsilon),...
        'M',num2str(pde.M),'S1=',num2str(pde.S1),'Nx=',num2str(Nx),'Ny=',num2str(Ny),'.txt'];
    T = table(dt_array,error);
    writetable(T,name);
end

%% results:
%% CAC_Vesicle_2D_MSAV_1st
% Accuracy test with exact solution.
%     dt     &   Error_L2   & Order
% 5.0000e-03 & 7.5023e+00 & 0.00
% 2.5000e-03 & 6.9938e+00 & 0.10
% 1.2500e-03 & 6.6741e+00 & 0.07
% 6.2500e-04 & 6.3468e+00 & 0.07
% 3.1250e-04 & 6.0005e+00 & 0.08
% 1.5625e-04 & 4.9419e+00 & 0.28
% 7.8125e-05 & 3.1117e+00 & 0.67
%
% Stabilized method
%     dt     &   Error_L2   & Order
% 5.0000e-03 & 3.7487e-02 & 0.00
% 2.5000e-03 & 1.6138e-02 & 1.22
% 1.2500e-03 & 7.4169e-03 & 1.12
% 6.2500e-04 & 3.5503e-03 & 1.06
% 3.1250e-04 & 1.7329e-03 & 1.03
% 1.5625e-04 & 8.5487e-04 & 1.02
% 7.8125e-05 & 4.2428e-04 & 1.01





