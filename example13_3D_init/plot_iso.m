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
para.C0 = 100; % SAV
para.epsilon = 6*pi/128;
% para.epsilon = 0.07;
para.M = 1;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

% Time: dt T
T = 5;
t0 = 0;
tsave = 0.04;

% dt_array = 0.01./2.^(1); 
dt_array = 0.01./2.^(2:9); 
% dt_array = 0.01./2.^(0:1);
dt_ref = 5e-4;

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
pdename = str2func('ex13_Vesicles_data_eight');
% pdename = str2func('ex13_Vesicles_data_six');
% pdename = str2func('ex13_Vesicles_data_five');
pdename = str2func('ex13_Vesicles_data_four');
% pdename = str2func('ex13_Vesicles_data_three_elliptic');
% pdename = str2func('ex13_Vesicles_data_threeOoo');
% pdename = str2func('ex13_Vesicles_data_twoOO');

scheme0 = 'CAC_Vesicle_3D_';

% scheme1_array = {'LM0_SAV_','LM1_SAV_','LM3_LM_'};
scheme1 = {'LM3_LM_'};

% scheme2_array = {'1st','BDF'};
scheme2 = {'BDF'};

dirname = [func2str(pdename),'_',[scheme1{1},scheme2{1}]];
datadir = [dirname,'/data'];
figdir  = [dirname,'/',dirname];

X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);
Z = load([datadir, '/Z.txt']);

figure(5)
% for t= [0]
% for t= [0:tsave:2] % 2.2:0.4:5]
% for t= [0 0.04 0.2 0.4 1.2 2]
for t= [0 0.12 0.8 2]
    t
    clf;
    ssp = [datadir '/phi_t=' num2str(t) '.txt'];
    phi = load(ssp);

%     subplot(1,2,1)
    X = reshape(X,Nx,Ny,Nz);
    Y = reshape(Y,Nx,Ny,Nz);
    Z = reshape(Z,Nx,Ny,Nz);
    phiu = reshape(phi,Nx,Ny,Nz);
    p0 = patch(isosurface(X,Y,Z,phiu,0));
    set(p0,'FaceColor','red', 'EdgeColor','none');

%     hold on
% 
%     filename = [datadir '/phi_b_t=' num2str(t)];
%     ss = [filename '.txt'];
%     phiu = load(ss);
%     phiu = reshape(phiu,Nx,Ny,Nz);    
%     p0 = patch(isosurface(xx,yy,zz,phiu,0.2));
%     set(p0,'FaceColor','blue', 'EdgeColor','none');

    daspect([1 1 1]); 
    camlight;  
    lighting phong;
%     axis image;
%     axis tight;
    axis square;
    axis equal;
    axis off;
%     axis([min(xx(:)) min(xx(:))min(yy(:)) max(yy(:)) min(zz(:)) max(zz(:)) ]);
%     axis 'auto xy';
    xlim([0, 2*pi])
    ylim([0, 2*pi])
    zlim([0, 2*pi])
    box on; 

    view(-161,16);
%     view(115,20)    
    
    ax = gca;
    c = ax.Color;
    ax.Color = [0.30,0.75,0.90];

%     xlabel('$x$','FontSize',18,'Interpreter','LaTex');
%     ylabel('$y$','FontSize',18,'Interpreter','LaTex');
%     zlabel('$z$','FontSize',18,'Interpreter','LaTex');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    drawnow;
    set(gcf, 'InvertHardCopy', 'off');
    set(0,'defaultfigurecolor','w') 
    figname = ['/Users/liq/work/05_AcademicActivity/2024_05_兰州大学/Vesicle_slide/figure_vesicle_01_Lagrange/' dirname '_phi_t=' num2str(t) '.png'];
    print(figname,'-dpng', '-r300')
end


