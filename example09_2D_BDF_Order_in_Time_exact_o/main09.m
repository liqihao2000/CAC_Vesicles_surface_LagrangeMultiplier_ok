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

N  = 128;
Nx = N;
Ny = N;

% Parameters
para.C0 = 1; % SAV
para.epsilon = 6*pi/128;
para.M = 1;
para.M1 = 50;
para.M2 = 50;

PDE = 'data2';

if 1 == strcmp(PDE,'data1')
    T = 0.0001;
    dt_array = 0.00001./2.^(0:4)';
    dt_ref = 1e-8;
    pde = ex03_Vesicles_data(para);
end

if 1 == strcmp(PDE,'data2')
    T = 0.02;
    dt_array = 0.001./2.^(0:6)';
    dt_ref = 1e-8;
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
    Vesicle_CH_2D_BDF_MSAV(pde,domain,Nx,Ny,time,option);
end
for k = 1:maxIt
    dt = dt_array(k);
    time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
    v2 = Vesicle_CH_2D_BDF_MSAV(pde,domain,Nx,Ny,time,option);
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
%    error(k,1) = sqrt(sum(sum((phi_exact - phi).^2))*hx*hy);   % L2
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

%% Vesicle_CH_2D_BDF_MSAV

%% ex03_Vesicles_data
% epsilon=0.147,t=0.0001/0.000, dt=6.25e-07, Nx=128, Ny=128, timeElapsed=3.433871
% dt     &   Error_L2	   &  Order 
% 1.0000e-05  4.665188e-10  0.00 
% 5.0000e-06  1.212498e-10  1.94 
% 2.5000e-06  3.086101e-11  1.97 
% 1.2500e-06  7.815288e-12  1.98 
% 6.2500e-07  2.014429e-12  1.96 

%% ex03_2_Vesicles_data
% epsilon=0.147,t=0.0200/0.020, dt=1.56e-05, Nx=128, Ny=128, timeElapsed=13.167627
% dt     &   Error_L2	   &  Order 
% 1.0000e-03  1.138941e-03  0.00 
% 5.0000e-04  1.828231e-04  2.64 
% 2.5000e-04  3.691805e-05  2.31 
% 1.2500e-04  8.333645e-06  2.15 
% 6.2500e-05  1.981074e-06  2.07 
% 3.1250e-05  4.829863e-07  2.04 
% 1.5625e-05  1.192404e-07  2.02 



