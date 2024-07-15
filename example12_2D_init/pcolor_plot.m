close all;
clear; clc;

% add path
addpath('../','-begin');

% Space: Domain and N
domain.left   =  0;
domain.right  =2*pi;
domain.bottom = -pi;
domain.top    =  pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

N = 64;
Nx = N;
Ny = N;

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
pdename = str2func('ex12_Vesicles_data_twoOO');
pdename = str2func('ex12_Vesicles_data_twoO_o');
% pdename = str2func('ex12_Vesicles_data_threeOoo');
% pdename = str2func('ex12_Vesicles_data_three_elliptic');
% pdename = str2func('ex12_Vesicles_data_fiveO');
% pdename = str2func('ex12_Vesicles_data_eightO');

scheme0 = 'CAC_Vesicle_2D_';

% scheme1_array = {'LM0_SAV_','LM1_SAV_','LM3_LM_'};
scheme1 = {'LM0_SAV_'};
scheme1 = {'LM3_LM_'};

% scheme2_array = {'1st','BDF'};
scheme2 = {'BDF'};

dirname = [func2str(pdename),'_',[scheme1{1},scheme2{1}]];
datadir = [dirname,'/data'];
figdir  = [dirname,'/',dirname];

X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);

figure(5)
% for t= [0]
for t= [0:tsave:2 2.2:0.4:5]
    t
    ssp = [datadir '/phi_t=' num2str(t) '.txt'];
    phi = load(ssp);

%     subplot(1,2,1)
    s=pcolor(X,Y,phi);
    s.FaceColor='interp';
    s.EdgeColor='interp';
%     view(0,90);
    colormap jet;
%     colormap viridis;
    axis square;
    axis tight;
    axis off;
%     subplot(1,2,2)
%     s=pcolor(lap_diff(phi,k2));
%     s.FaceColor='interp';
%     s.EdgeColor='interp';
% %     view(0,90);
%     colormap gray;
%     axis square;
%     axis tight;
%     axis off;
%     colorbar('Position',[0.845 0.18 0.03 0.66],'Fontsize',15);
%     caxis([-0.7 0.7])
%     colorbar off;
    
    drawnow;
    
%     figname = [figdir '/phi_t=' num2str(t) '.eps'];
%     print(figname,'-depsc2', '-r120')

    figname = [figdir '_phi_t=' num2str(t) '.png'];
%     print(figname,'-dpng', '-r300')
end


