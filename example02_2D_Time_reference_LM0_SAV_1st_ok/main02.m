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


%% Run:
delete *.mat
if ~isfield(pde,'exact') || ~isfield(pde,'rhs')
    time = struct('T',T,'t0',t0,'dt',dt_ref,'tsave',tsave);
    CAC_Vesicle_2D_LM0_SAV_1st(pde,domain,Nx,Ny,time,option);
%     CAC_Vesicle_2D_LM1_SAV_1st(pde,domain,Nx,Ny,time,option);
%     CAC_Vesicle_2D_LM3_LM_1st(pde,domain,Nx,Ny,time,option);
end
for k = 1:maxIt
    dt = dt_array(k);
    time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
    v2 = CAC_Vesicle_2D_LM0_SAV_1st(pde,domain,Nx,Ny,time,option);
%     v2 = CAC_Vesicle_2D_LM1_SAV_1st(pde,domain,Nx,Ny,time,option);
%     v2 = CAC_Vesicle_2D_LM3_LM_1st(pde,domain,Nx,Ny,time,option);
end

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
else
    Lx = domain.right - domain.left;
    Ly = domain.top   - domain.bottom;
    hx = Lx/Nx;
    hy = Ly/Ny;
    x  = domain.left   + hx*(0:Nx-1);
    y  = domain.bottom + hy*(0:Ny-1);
    [xx,yy] = ndgrid(x,y);
    phi_exact = pde.exact(xx,yy,T);
end
for k = 1:maxIt
    dt = dt_array(k);
    name=['phi_e',num2str(pde.epsilon),'M',num2str(pde.M),...
        'Nx=',num2str(Nx),'Ny=',num2str(Ny),'dt=',num2str(dt)];
    filenamek=[name '.mat'];
    load(filenamek,'phi','hx','hy');
    err    = fft2((phi_exact - phi).^2);
    error(k,1) = sqrt(err(1,1)*hx*hy);   % L2
    clear phi;
end
order(2:maxIt) = log(error(1:maxIt-1)./error(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));

%% Display error and order
fprintf('dt     &   Error_L2\t   &  Order \n');
for k = 1:maxIt
    fprintf('%.4e  %.6e  %.2f \n',dt_array(k),error(k),order(k));
    %         fprintf('1/%d\t& %.6e\t& %.4f %s \n',1./dt_array(k),error(k),order(k),'\\');
end
fprintf('\n')

%% Plot
figure(2)
hh=loglog(dt_array,error,'*-');
xlabel('Time step $\delta t$','Interpreter','latex');
%  ylabel('Error\_max')
ylabel('$L^2$ error','Interpreter','latex');
grid on;
hold on;

%% Save error and order
name=['phi_e',num2str(para.epsilon),...
      'M',num2str(pde.M),'Nx=',num2str(N),'Ny=',num2str(N),'.txt'];
T = table(dt_array,error);
writetable(T,name);

%% results:
% M = 2;
%% CAC_Vesicle_2D_LM0_SAV_BDF
% ex03_1_Vesicles_data
% lambda=-1.8754e+02,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=2.186915
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  1.220740e-03  0.00 
% 5.0000e-06  2.883158e-04  2.08 
% 2.5000e-06  7.022521e-05  2.04 
% 1.2500e-06  1.733285e-05  2.02 
% 6.2500e-07  4.304832e-06  2.01 
% 3.1250e-07  1.071948e-06  2.01 
% 1.5625e-07  2.667349e-07  2.01 

% ex03_2_Vesicles_data
% lambda=-2.3895e+01,epsilon=0.147,t=0.00010/0.0001, dt=1.53e-09, Nx=32, Ny=32, timeElapsed=112.782871
% dt     &   Error_L2	   &  Order 
% 3.1250e-06  2.073358e-04  0.00 
% 1.5625e-06  5.414485e-05  1.94 
% 7.8125e-07  1.384476e-05  1.97 
% 3.9063e-07  3.501126e-06  1.98 
% 1.9531e-07  8.803511e-07  1.99 
% 9.7656e-08  2.207121e-07  2.00 
% 4.8828e-08  5.524105e-08  2.00 
% 2.4414e-08  1.380284e-08  2.00 
% 1.2207e-08  3.434740e-09  2.01 
% 6.1035e-09  8.425099e-10  2.03 
% 3.0518e-09  1.972388e-10  2.09 
% 1.5259e-09  4.464986e-11  2.14 

%% CAC_Vesicle_2D_LM1_SAV_BDF
% ex03_1_Vesicles_data
% lambda=-1.8754e+02,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=2.676092
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  1.220740e-03  0.00 
% 5.0000e-06  2.883158e-04  2.08 
% 2.5000e-06  7.022521e-05  2.04 
% 1.2500e-06  1.733285e-05  2.02 
% 6.2500e-07  4.304832e-06  2.01 
% 3.1250e-07  1.071948e-06  2.01 
% 1.5625e-07  2.667349e-07  2.01 

% ex03_2_Vesicles_data
% lambda=-2.3895e+01,epsilon=0.147,t=0.00010/0.0001, dt=1.53e-09, Nx=32, Ny=32, timeElapsed=135.433879
% dt     &   Error_L2	   &  Order 
% 3.1250e-06  2.073358e-04  0.00 
% 1.5625e-06  5.414485e-05  1.94 
% 7.8125e-07  1.384476e-05  1.97 
% 3.9063e-07  3.501126e-06  1.98 
% 1.9531e-07  8.803511e-07  1.99 
% 9.7656e-08  2.207121e-07  2.00 
% 4.8828e-08  5.524105e-08  2.00 
% 2.4414e-08  1.380284e-08  2.00 
% 1.2207e-08  3.434740e-09  2.01 
% 6.1035e-09  8.425099e-10  2.03 
% 3.0518e-09  1.972388e-10  2.09 
% 1.5259e-09  4.464986e-11  2.14 

%% CAC_Vesicle_2D_LM3_LM_BDF
% ex03_1_Vesicles_data
% lambda=2.5410e+01,eta=1.0000e+00,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=8.517206
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  2.071581e-05  0.00 
% 5.0000e-06  5.193254e-06  2.00 
% 2.5000e-06  1.386889e-06  1.90 
% 1.2500e-06  3.474933e-07  2.00 
% 6.2500e-07  1.754235e-07  0.99 
% 3.1250e-07  4.384467e-08  2.00 
% 1.5625e-07  1.093145e-08  2.00 

% ex03_2_Vesicles_data
% lambda=-2.4281e+01,eta=1.0000e+00,epsilon=0.147,t=0.00010/0.0001, dt=1.53e-09, Nx=32, Ny=32, timeElapsed=452.025976
% dt     &   Error_L2	   &  Order 
% 3.1250e-06  2.685011e-04  0.00 
% 1.5625e-06  2.129919e-04  0.33 
% 7.8125e-07  2.176701e-04  -0.03 
% 3.9063e-07  2.205794e-04  -0.02 
% 1.9531e-07  2.227479e-04  -0.01 
% 9.7656e-08  2.313768e-04  -0.05 
% 4.8828e-08  6.305810e-04  -1.45 
% 2.4414e-08  5.551141e-04  0.18 
% 1.2207e-08  3.528401e-04  0.65 
% 6.1035e-09  5.445027e-05  2.70 
% 3.0518e-09  1.921924e-10  18.11 
% 1.5259e-09  3.876093e-11  2.31 




%% exp data
% lambda=-3.3318e+01,eta=1.0000e+00,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=7.760688
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  1.159633e-06  0.00 
% 5.0000e-06  2.898893e-07  2.00 
% 2.5000e-06  7.252109e-08  2.00 
% 1.2500e-06  1.813498e-08  2.00 
% 6.2500e-07  4.532935e-09  2.00 
% 3.1250e-07  1.131745e-09  2.00 
% 1.5625e-07  2.813933e-10  2.01 


