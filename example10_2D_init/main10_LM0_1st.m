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
para.M  = 1;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

para.S1 = 0;
para.S2 = 0;
para.S3 = 0;

para.name = 'ex10_data_LM0_1st';

% Time: dt T
T = 20;
t0 = 0;
tsave = 0.02*T;

dt_array = 0.01./2.^(0:4);
dt_ref = 4e-4;

maxIt = length(dt_array);

%% option
option.plotflag  = 0;
option.printflag = 1;
option.vtkflag  = 0;
option.saveflag  = 0;
option.savefinal  = 0;
option.energyflag = 1;
option.tol = 1e-14;
option.tolit = 1e-11;
option.maxit = 2000;

pde = ex10_Vesicles_data(para);

% %% Run:
% for k = 1:maxIt
%     dt = dt_array(k);
%     time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
%     CAC_Vesicle_2D_LM0_SAV_1st(pde,domain,Nx,Ny,time,option);
% %     CAC_Vesicle_2D_LM1_SAV_1st(pde,domain,Nx,Ny,time,option);
% %     CAC_Vesicle_2D_LM3_LM_1st(pde,domain,Nx,Ny,time,option);
% end

%% plot
% energy
figure(1);
hold on;
% lineType = {'.-', 's-','*-','o-','+-' ,'.-' ,'--k',':','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};

lineType = {'-','-','-','-' ,'-','b-','-','-','-','-k','-b'};
for k = 1:maxIt
    energy=load([pde.name,'_dt_',num2str(dt_array(k)),'_energy.txt']);
    tmp = 1;
    plot(energy(tmp:1:end,1),energy(tmp:1:end,3),char(lineType(k)),'LineWidth',2.5);
end
h = legend('$1:\delta t = 0.01/2^2$','$2:\delta t = 0.01/2^3$','$3:\delta t = 0.01/2^4$',...
    '$4:\delta t = 0.01/2^5$','$5:\delta t = 0.01/2^6$','$6:\delta t = 0.01/2^7$');
% set(h,'interpreter','latex');
xlabel('Time','Fontsize',24);ylabel('Energy $R$','Fontsize',24,'interpreter','latex');
% xlim([0,5])
% % ylim([6490,6650])
% % yticks([0:2:20])
set(h,'box','off','interpreter','latex');
set(gca,'FontSize',22);
set(gca,'linewidth',1.8)
% set(gca,'xtick',0:0.2:T);
% set(gca,'ytick',3:0.1:4);
grid on;
box on;
figure_FontSize=35;
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰

% figname1 = ['../../../../../papers/paper15_Vesicles/01_CHNS_cn/1Revision/figure_Vesicles_eps/',pde.name,'_stability_0','.eps']
% print(figname1,'-depsc', '-r250')

figname1 = ['/Users/liq/work/10_UndergraduateThesis/2024/2024_chenhuiyi/CHD_Bachelor_ChenHuiyi/figs/',pde.name,'_stability_0','.png'];
% print(figname1,'-dpng', '-r250')

% mass
figure(2);
hold on;
% lineType = {'.-', 's-','*-','o-','+-' ,'.-' ,'--k',':','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};
% lineType = {'>-', 's-','*-','o-','+-' ,'.-' ,'--k','-b'};

lineType = {'-','-','-','-' ,'-','b-','-','-','-','-k','-b'};
for k = 1:maxIt
    mass=load([pde.name,'_dt_',num2str(dt_array(k)),'_mass.txt']);
%     energy=load([pde.name,'_dt_',num2str(dt_array(k)),'_mass.txt']);
    tmp = 1;
    plot(mass(tmp:1:end,1),mass(tmp:1:end,2),char(lineType(k)),'LineWidth',2.5);
end
h = legend('$1:\delta t = 0.01/2^2$','$2:\delta t = 0.01/2^3$','$3:\delta t = 0.01/2^4$',...
    '$4:\delta t = 0.01/2^5$','$5:\delta t = 0.01/2^6$','$6:\delta t = 0.01/2^7$');
% set(h,'interpreter','latex');
xlabel('Time','Fontsize',24);ylabel('Surface','Fontsize',24,'interpreter','latex');
% xlim([0,5])
% % ylim([6490,6650]
% % yticks([0:2:20])
set(h,'box','off','Location','southeast','interpreter','latex');
set(gca,'FontSize',22);
set(gca,'linewidth',1.8)
% set(gca,'xtick',0:0.2:T);
% set(gca,'ytick',3:0.1:4);
grid on;
box on;

figure_FontSize=35;
set(get(gca,'XLabel'),'FontSize',figure_FontSize);
set(get(gca,'YLabel'),'FontSize',figure_FontSize);
set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰

% figname2 = ['../../../../../papers/paper15_Vesicles/01_CHNS_cn/1Revision/figure_Vesicles_eps/',pde.name,'_surface_0','.eps'];
% print(figname2,'-depsc', '-r250')

figname1 = ['/Users/liq/work/10_UndergraduateThesis/2024/2024_chenhuiyi/CHD_Bachelor_ChenHuiyi/figs/',pde.name,'_surface_0','.png'];
% print(figname1,'-dpng', '-r250')



