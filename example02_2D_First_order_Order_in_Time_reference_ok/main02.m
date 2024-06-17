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
    T = 0.0002;
    dt_array = 0.00001./2.^(0:6)';
    dt_ref = 1e-8;
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
%     CAC_Vesicle_2D_LM0_SAV_1st(pde,domain,Nx,Ny,time,option);
%     CAC_Vesicle_2D_LM1_SAV_1st(pde,domain,Nx,Ny,time,option);
    CAC_Vesicle_2D_LM3_LM_1st(pde,domain,Nx,Ny,time,option);
end
for k = 1:maxIt
    dt = dt_array(k);
    time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
%     v2 = CAC_Vesicle_2D_LM0_SAV_1st(pde,domain,Nx,Ny,time,option);
%     v2 = CAC_Vesicle_2D_LM1_SAV_1st(pde,domain,Nx,Ny,time,option);
    v2 = CAC_Vesicle_2D_LM3_LM_1st(pde,domain,Nx,Ny,time,option);
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
name=['phi_e',num2str(para.epsilon),'M',num2str(para.M),...
      'Nx=',num2str(N),'Ny=',num2str(N)];
fileID = fopen([name,'.txt'],'w');
% fprintf(fileID,'%6s\n','%% Results');
% fprintf(fileID,'%6s\n','% dt	   &   Error_L2	   &  Order');
% A = [dt_array error];
% fprintf(fileID,'%.12f   %.4e   \n',A');
fprintf(fileID,'%.12f     %.4e      %.2f \n',[dt_array,error,order]');
fclose(fileID);

%% results:
%% M = 2;
%% CAC_Vesicle_2D_LM0_SAV_1st
% ex03_1_Vesicles_data
% lambda=-1.8762e+02,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=2.394364
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  8.143822e-03  0.00 
% 5.0000e-06  4.012460e-03  1.02 
% 2.5000e-06  1.988420e-03  1.01 
% 1.2500e-06  9.867918e-04  1.01 
% 6.2500e-07  4.885687e-04  1.01 
% 3.1250e-07  2.401038e-04  1.02 
% 1.5625e-07  1.160329e-04  1.05 

% ex03_2_Vesicles_data
% lambda=-2.3894e+01,epsilon=0.147,t=0.00010/0.0001, dt=1.53e-09, Nx=32, Ny=32, timeElapsed=126.373443
% dt     &   Error_L2	   &  Order 
% 3.1250e-06  1.164541e-03  0.00 
% 1.5625e-06  5.861618e-04  0.99 
% 7.8125e-07  2.939122e-04  1.00 
% 3.9063e-07  1.470221e-04  1.00 
% 1.9531e-07  7.338576e-05  1.00 
% 9.7656e-08  3.651981e-05  1.01 
% 4.8828e-08  1.807489e-05  1.01 
% 2.4414e-08  8.849439e-06  1.03 
% 1.2207e-08  4.235970e-06  1.06 
% 6.1035e-09  1.929049e-06  1.13 
% 3.0518e-09  7.755418e-07  1.31 
% 1.5259e-09  1.987766e-07  1.96 

%% CAC_Vesicle_2D_LM1_SAV_1st
% ex03_1_Vesicles_data
% lambda=-1.8762e+02,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=2.517449
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  5.059063e-03  0.00 
% 5.0000e-06  4.012460e-03  0.33 
% 2.5000e-06  1.988420e-03  1.01 
% 1.2500e-06  9.867918e-04  1.01 
% 6.2500e-07  4.885687e-04  1.01 
% 3.1250e-07  2.401038e-04  1.02 
% 1.5625e-07  1.160329e-04  1.05 

% ex03_2_Vesicles_data
% lambda=-2.3894e+01,epsilon=0.147,t=0.00010/0.0001, dt=1.53e-09, Nx=32, Ny=32, timeElapsed=191.899195
% dt     &   Error_L2	   &  Order 
% 3.1250e-06  1.164541e-03  0.00 
% 1.5625e-06  5.861618e-04  0.99 
% 7.8125e-07  2.939122e-04  1.00 
% 3.9063e-07  1.470221e-04  1.00 
% 1.9531e-07  7.338576e-05  1.00 
% 9.7656e-08  3.651981e-05  1.01 
% 4.8828e-08  1.807489e-05  1.01 
% 2.4414e-08  8.849439e-06  1.03 
% 1.2207e-08  4.235970e-06  1.06 
% 6.1035e-09  1.929049e-06  1.13 
% 3.0518e-09  7.755418e-07  1.31 
% 1.5259e-09  1.987766e-07  1.96 

%% CAC_Vesicle_2D_LM3_LM_1st
% ex03_1_Vesicles_data
% lambda=2.5432e+01,eta=1.0000e+00,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=4.738573
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  2.868115e-04  0.00 
% 5.0000e-06  1.431607e-04  1.00 
% 2.5000e-06  7.141278e-05  1.00 
% 1.2500e-06  3.550577e-05  1.01 
% 6.2500e-07  2.977461e-05  0.25 
% 3.1250e-07  1.464552e-05  1.02 
% 1.5625e-07  7.080751e-06  1.05 

% ex03_2_Vesicles_data
% lambda=-2.4281e+01,eta=1.0000e+00,epsilon=0.147,t=0.00010/0.0001, dt=1.53e-09, Nx=32, Ny=32, timeElapsed=382.216807
% dt     &   Error_L2	   &  Order 
% 3.1250e-06  1.133841e-03  0.00 
% 1.5625e-06  5.718248e-04  0.99 
% 7.8125e-07  3.188446e-04  0.84 
% 3.9063e-07  2.313493e-04  0.46 
% 1.9531e-07  2.199423e-04  0.07 
% 9.7656e-08  6.179738e-04  -1.49 
% 4.8828e-08  5.482338e-04  0.17 
% 2.4414e-08  3.491430e-04  0.65 
% 1.2207e-08  5.296901e-05  2.72 
% 6.1035e-09  1.902867e-06  4.80 
% 3.0518e-09  7.650172e-07  1.31 
% 1.5259e-09  1.960811e-07  1.96 


%% exp data
% lambda=-3.3318e+01,eta=1.0000e+00,epsilon=0.147,t=0.00020/0.0002, dt=1.56e-07, Nx=32, Ny=32, timeElapsed=4.997044
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  1.529054e-05  0.00 
% 5.0000e-06  7.647791e-06  1.00 
% 2.5000e-06  3.821975e-06  1.00 
% 1.2500e-06  1.904747e-06  1.00 
% 6.2500e-07  9.450486e-07  1.01 
% 3.1250e-07  4.649283e-07  1.02 
% 1.5625e-07  2.248003e-07  1.05 



