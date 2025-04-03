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

N  = 32;
Nx = N;
Ny = N;
Nz = N;

% Parameters
para.epsilon = 6*pi/128;
para.M = 2;

% delete *.mat
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

    PDE = 'data1';

    if 1 == strcmp(PDE,'data1')
        T = 0.002;
        dt_array = 0.0001./2.^(0:6)';
        dt_ref = 1e-7;
        para.C0 = 100; % SAV
        pde = ex02_3D_Vesicles_data(para);
    end

    % Time: dt T

    t0 = 0;
    tsave = 0.2*T;

    maxIt = length(dt_array);

    %% option
    option.plotflag  = 0;
    option.printflag = 1;
    option.vtkflag  = 0;
    option.saveflag  = 0;
    option.savefinal  = 1;
    option.energyflag = 0;
    option.tol = 1e-14;
    option.tolit = 1e-11;
    option.maxit = 2000;

%     %% Run:
%     %     delete *.mat
%     if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
%         time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
% %         CAC_Vesicle_3D_LM0_SAV_BDF(pde,domain,Nx,Ny,Nz,time,option);
%         CAC_Vesicle_3D_LM1_SAV_BDF(pde,domain,Nx,Ny,Nz,time,option);
% %         CAC_Vesicle_3D_LM3_LM_BDF(pde,domain,Nx,Ny,Nz,time,option);
%     end
%     for k = 1:maxIt
%         dt = dt_array(k);
%         time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
% %         v2 = CAC_Vesicle_3D_LM0_SAV_BDF(pde,domain,Nx,Ny,Nz,time,option);
%         v2 = CAC_Vesicle_3D_LM1_SAV_BDF(pde,domain,Nx,Ny,Nz,time,option);
% %         v2 = CAC_Vesicle_3D_LM3_LM_BDF(pde,domain,Nx,Ny,Nz,time,option);
%     end

    %% Compute order of convergence
    error=zeros(maxIt,1);
    order=zeros(maxIt,1);
    if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
        name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),'S1=',num2str(pde.S1),...
            'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt_ref)];
        filename=[name '.mat'];
        load(filename,'phi');
        phi_exact = phi;
        clear phi;
        fprintf('Accuracy test with reference solution.\n');
    end
    for k = 1:maxIt
        dt = dt_array(k);
        name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),'S1=',num2str(pde.S1),...
            'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'dt=',num2str(dt)];
        filenamek=[name '.mat'];
        load(filenamek,'phi','xx','yy','zz','hx','hy','hz');

        if ~exist('phi_exact','var') && ~exist('psi_exact','var')
            phi_exact = pde.exactphi(xx,yy,zz,T);
            fprintf('Accuracy test with exact solution.\n');
        end

        err    = fftn((phi_exact - phi).^2);
        error(k,1) = sqrt(err(1,1,1)*hx*hy);   % L2
        clear phi;
    end
    order(2:maxIt) = log(error(1:maxIt-1)./error(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

    %% Display error and order
    if 0 ~= para.S1
        fprintf('Stabilized method\n');
    end
    fprintf('    dt     &   Error_L2 & Order \n');
    for k = 1:maxIt
        fprintf('%.4e & %.4e & %.2f \n',dt_array(k),error(k),order(k));
        %         fprintf('1/%d\t& %.6e\t& %.4f %s \n',1./dt_array(k),error(k),order(k),'\\');
    end
    fprintf('\n')

    % %% Plot
    % figure(2)
    % hh=loglog(dt_array,error,'*-');
    % xlabel('Time step $\delta t$','Interpreter','latex');
    % %  ylabel('Error\_max')
    % ylabel('$L^2$ error','Interpreter','latex');
    % grid on;
    % hold on;

    %% Save error and order
    name=['phi_e',num2str(para.epsilon),...
        'M',num2str(pde.M),'S1=',num2str(pde.S1),'Nx=',num2str(Nx),'Ny=',num2str(Ny),'Nz=',num2str(Nz),'.txt'];
    T = table(dt_array,error);
    writetable(T,name);
end

%% results:
% lambda=-3.3197e+01,epsilon=0.147,t=0.00200/0.0020, dt=1.56e-06, Nx=32, Ny=32, Nz=32, timeElapsed=73.719963
% Accuracy test with reference solution.
%     dt     &   Error_L2 & Order 
% 1.0000e-04 & 5.3078e-05 & 0.00 
% 5.0000e-05 & 1.3416e-05 & 1.98 
% 2.5000e-05 & 3.3726e-06 & 1.99 
% 1.2500e-05 & 8.4548e-07 & 2.00 
% 6.2500e-06 & 2.1163e-07 & 2.00 
% 3.1250e-06 & 5.2901e-08 & 2.00 
% 1.5625e-06 & 1.3187e-08 & 2.00 

% Accuracy test with reference solution.
% Stabilized method
%     dt     &   Error_L2 & Order 
% 1.0000e-04 & 4.8888e-04 & 0.00 
% 5.0000e-05 & 1.2944e-04 & 1.92 
% 2.5000e-05 & 3.3559e-05 & 1.95 
% 1.2500e-05 & 8.5676e-06 & 1.97 
% 6.2500e-06 & 2.1660e-06 & 1.98 
% 3.1250e-06 & 5.4431e-07 & 1.99 
% 1.5625e-06 & 1.3607e-07 & 2.00 






