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

N  = 32;
Nx = N;
Ny = N;

% Parameters
para.epsilon = 6*pi/128;
para.M = 2;

PDE = 'data1';

if 1 == strcmp(PDE,'data1')
    T = 0.002;
    dt_array = 0.0001./2.^(0:6)';
    dt_ref = 1e-7;
    para.C0 = 100; % SAV
    pde = ex03_1_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.0001;
    dt_array = 0.0001./2.^(5:16)';
    dt_ref = 1e-9;
    para.C0 = 1000; % SAV
    pde = ex03_2_Vesicles_data(para);
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

% %% Run:
% delete *.mat
% if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
%     time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
%     CAC_Vesicle_2D_LM0_SAV_BDF(pde,domain,Nx,Ny,time,option);
% %     CAC_Vesicle_2D_LM1_SAV_BDF(pde,domain,Nx,Ny,time,option);
% %     CAC_Vesicle_2D_LM3_LM_BDF(pde,domain,Nx,Ny,time,option);
% end
% for k = 1:maxIt
%     dt = dt_array(k);
%     time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
%     v2 = CAC_Vesicle_2D_LM0_SAV_BDF(pde,domain,Nx,Ny,time,option);
% %     v2 = CAC_Vesicle_2D_LM1_SAV_BDF(pde,domain,Nx,Ny,time,option);
% %     v2 = CAC_Vesicle_2D_LM3_LM_BDF(pde,domain,Nx,Ny,time,option);
% end

%% Compute order of convergence
error=zeros(maxIt,1);
order=zeros(maxIt,1);
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
          'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt_ref)];
    filename=[name '.mat'];
    load(filename,'phi');
    phi_exact = phi;
    clear phi;
    fprintf('Accuracy test with reference solution.\n');
end
for k = 1:maxIt
    dt = dt_array(k);
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filenamek=[name '.mat'];
    load(filenamek,'phi','xx','yy','hx','hy');

    if ~exist('phi_exact','var') && ~exist('psi_exact','var') 
            phi_exact = pde.exactphi(xx,yy,T);
            fprintf('Accuracy test with exact solution.\n');
    end

    err    = fft2((phi_exact - phi).^2);
    error(k,1) = sqrt(err(1,1)*hx*hy);   % L2
    clear phi;
end
order(2:maxIt) = log(error(1:maxIt-1)./error(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

%% Display error and order
fprintf('    dt     &   Error_L2   & Order \n');
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
      'M',num2str(pde.M),'Nx=',num2str(N),'Ny=',num2str(N),'.txt'];
T = table(dt_array,error);
writetable(T,name);

%% results:
% Accuracy test with reference solution.
%     dt     &   Error_L2   & Order 
% 1.0000e-04 & 1.4625e-04 & 0.00 
% 5.0000e-05 & 3.7117e-05 & 1.98 
% 2.5000e-05 & 9.3477e-06 & 1.99 
% 1.2500e-05 & 2.3454e-06 & 1.99 
% 6.2500e-06 & 5.8729e-07 & 2.00 
% 3.1250e-06 & 1.4684e-07 & 2.00 
% 1.5625e-06 & 3.6613e-08 & 2.00 



