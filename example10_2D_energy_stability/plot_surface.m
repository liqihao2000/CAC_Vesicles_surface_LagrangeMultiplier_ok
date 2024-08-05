clear; clc;
close all;

% Parameters
para.C0 = 100; % SAV
para.epsilon = 6*pi/128;
% para.epsilon = 0.07;
para.M = 1;
para.S1 = 4;
para.S2 = 4;
para.S3 = 1;

% Time: dt T
T = 2;
t0 = 0;
tsave = 0.5*T;

% dt_array = 0.01./2.^(3:9);
dt_array = 0.01./2.^(4);
dt_ref = 4e-4;

maxIt = length(dt_array);

%% user parameters
pdename = str2func('ex10_Vesicles_data');
scheme0 = 'CAC_Vesicle_2D_';

scheme1_array = {'MSAV_','LM0_SAV_','LM1_SAV_','LM3_LM_'};
scheme1_array = {'MSAV_','LM1_SAV_','LM3_LM_'};

% scheme2_array = {'1st','BDFs'};
scheme2_array = {'BDF'};

legend_str =string(['$\delta t = ',num2str(dt_array(1)),'$']);
for kk_dt = 2:length(dt_array)
    legend_str= [legend_str,['$\delta t = ',num2str(dt_array(kk_dt)),'$']];
end

%% surface Area
% index_fig = 1;
figure(1)
for i_2 = 1:length(scheme2_array)
    for i_1 = 1:length(scheme1_array)

        %         figure(index_fig);
        %         index_fig = index_fig + 1;

        scheme1 = scheme1_array{i_1};
        scheme2 = scheme2_array{i_2};

        dirname = [func2str(pdename),'_',[scheme1,scheme2]];

        scheme = [scheme0,scheme1,scheme2];
        para.name = [func2str(pdename),[scheme1,scheme2]];

        pde = ex10_Vesicles_data(para);
        hold on;
        for k = 1:maxIt
            figname_mass = [dirname,'_S1_',num2str(pde.S1),'_dt_',num2str(dt_array(k)),'_mass.txt'];

            mass=load(figname_mass);
            tmp = 1;
            plot(mass(tmp:1:end,1),abs((mass(tmp:1:end,3)-mass(1,3))./mass(1,3)),'-','LineWidth',2.5);
%             plot(mass(tmp:1:end,1),mass(tmp:1:end,3),'-','LineWidth',2.5);
        end
    end
end

% h = legend('MSAV','LM0SAV','LM1SAV','LM3LM');
h = legend('MSAV','SAV-BDF','LM-BDF');
xlabel('Time','Fontsize',18,'interpreter','latex');
ylabel('Ratio of Surface Area Difference','Fontsize',18,'interpreter','latex');
set(gca,'FontSize',24);
% set(gca,'linewidth',1.5);
% % xlim([0 18])
%         ylim([2 16])
set(h,'box','off','interpreter','latex','FontSize',16);
set(gca,'XTick',0:0.5:2)
% set(gca,'YTick',-25:5:10)
% % annotation('arrow',[0.435,0.35],[0.44,0.39],'LineWidth',2)
box on;
grid on;
% % set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.32]);
% % set (gca,'position',[0.15,0.11,0.8,0.8] );

set(gcf, 'InvertHardCopy', 'off');
set(0,'defaultfigurecolor','w') 
figname = ['/Users/liq/work/05_AcademicActivity/2024_05_兰州大学/Vesicle_slide/figure_vesicle_01_Lagrange/' dirname 'surface_ratio' '.png'];
print(figname,'-dpng', '-r300')


