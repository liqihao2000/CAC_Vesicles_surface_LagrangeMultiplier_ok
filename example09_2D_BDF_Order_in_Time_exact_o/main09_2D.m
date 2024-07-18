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
        para.S1 = 10;
        para.S2 = 10;
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
        T = 0.01;
        %     dt_array = 0.001./2.^(0:6)';
        dt_array = 0.0005./2.^(0:9)';
        dt_ref = 1e-8;
        para.C0 = 100; % SAV
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
        CAC_Vesicle_2D_MSAV_BDF(pde,domain,Nx,Ny,time,option);
        %     CAC_Vesicle_2D_LM0_SAV_BDF(pde,domain,Nx,Ny,time,option);
        %     CAC_Vesicle_2D_LM1_SAV_BDF(pde,domain,Nx,Ny,time,option);
        %     CAC_Vesicle_2D_LM3_LM_BDF(pde,domain,Nx,Ny,time,option);
    end
    for k = 1:maxIt
        dt = dt_array(k);
        time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
        v2 = CAC_Vesicle_2D_MSAV_BDF(pde,domain,Nx,Ny,time,option);
        %     v2 = CAC_Vesicle_2D_LM0_SAV_BDF(pde,domain,Nx,Ny,time,option);
        %     v2 = CAC_Vesicle_2D_LM1_SAV_BDF(pde,domain,Nx,Ny,time,option);
        %     v2 = CAC_Vesicle_2D_LM3_LM_BDF(pde,domain,Nx,Ny,time,option);
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
%% CAC_Vesicle_2D_MSAV_BDF
% Accuracy test with exact solution.
%     dt     &   Error_L2   & Order 
% 5.0000e-04 & 7.3919e-02 & 0.00 
% 2.5000e-04 & 7.7414e-02 & -0.07 
% 1.2500e-04 & 8.2347e-02 & -0.09 
% 6.2500e-05 & 8.5048e-02 & -0.05 
% 3.1250e-05 & 8.3158e-02 & 0.03 
% 1.5625e-05 & 8.2286e-02 & 0.02 
% 7.8125e-06 & 8.0692e-02 & 0.03 
% 3.9063e-06 & 7.7229e-02 & 0.06 
% 1.9531e-06 & 3.0787e-02 & 1.33 
% 9.7656e-07 & 1.3069e-11 & 31.13 
% 
% Stabilized method
%     dt     &   Error_L2   & Order 
% 5.0000e-04 & 2.3666e-04 & 0.00 
% 2.5000e-04 & 3.3662e-08 & 12.78 
% 1.2500e-04 & 9.6320e-07 & -4.84 
% 6.2500e-05 & 2.9723e-07 & 1.70 
% 3.1250e-05 & 7.8841e-08 & 1.91 
% 1.5625e-05 & 2.0146e-08 & 1.97 
% 7.8125e-06 & 5.0834e-09 & 1.99 
% 3.9063e-06 & 1.2762e-09 & 1.99 
% 1.9531e-06 & 3.1971e-10 & 2.00 
% 9.7656e-07 & 8.0004e-11 & 2.00 





