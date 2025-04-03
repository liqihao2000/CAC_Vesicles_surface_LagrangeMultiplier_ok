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

% dt_array = 0.01./2.^(1:7);
dt_array = 0.01./2.^(1:6);
dt_ref = 4e-4;

maxIt = length(dt_array);

%% user parameters
pdename = 'ex13_Vesicles_data_eight_';
scheme0 = 'CAC_Vesicle_3D_';

%     scheme1_array = {'MSAV_'};
scheme1_array = {'MSAV_','LM0_SAV_','LM1_SAV_','LM3_LM_','LM4_EX_LM_'};
scheme1_array = {'MSAV_'};
scheme1_array = {'LM4_EX_LM_'};

scheme1_array = {'MSAV_','LM4_EX_LM_'};


% scheme2_array = {'1st','BDF'};
scheme2_array = {'BDF'};

% dirname = [pdename,'_',[scheme1_array,scheme2_array]];
% datadir = [dirname,'/data'];
% figdir  = [dirname,'/',dirname];

pow = -log(dt_array*100)/log(2);
legend_str = string(['$\delta t = 0.01/2^{',num2str(pow(1)),'}$']);
for kk_dt = 2:length(dt_array)
    %     legend_str= [legend_str,['$\delta t = 0.01/2^{',num2str(dt_array(kk_dt)),'}$']];
    legend_str= [legend_str,['$\delta t = 0.01/2^{',num2str(pow(kk_dt)),'}$']];
end

figure_FontSize = 30;

index_fig = 1;


% energy
for i_2 = 1:length(scheme2_array)
    for i_1 = 1:length(scheme1_array)

        figure(index_fig);
        index_fig = index_fig + 1;

        scheme1 = scheme1_array{i_1};
        scheme2 = scheme2_array{i_2};

        scheme = [scheme0,scheme1,scheme2];
        para.name = [pdename,[scheme1,scheme2]];

        pde = ex13_Vesicles_data_eight(para);
        hold on;
        for k = 1:maxIt
            figname_energy = [pde.name,'_S1_',num2str(pde.S1),'_dt_',num2str(dt_array(k)),'_energy.txt'];

            energy=load(figname_energy);
            tmp = 1;
            plot(energy(tmp:1:end,1),energy(tmp:1:end,2),'-','LineWidth',2.5);
        end
        h = legend(legend_str);
        xlabel('Time','Fontsize',figure_FontSize,'interpreter','latex');
        ylabel('Discrete Energy','Fontsize',figure_FontSize-6,'interpreter','latex');
        set(gca,'FontSize',figure_FontSize-6);
        % set(gca,'linewidth',1.5);
        % % xlim([0 18])
        %         ylim([2 16])
        set(h,'box','off','interpreter','latex','FontSize',figure_FontSize-8);
        %         set(gca,'XTick',0:0.5:2)
        % set(gca,'YTick',-25:5:10)
        % % annotation('arrow',[0.435,0.35],[0.44,0.39],'LineWidth',2)
        box on;
        grid on;

        % 获取当前坐标轴
        ax = gca;

        % 设置坐标轴的边框粗细
        ax.Box = 'on'; % 确保边框打开
        ax.LineWidth = 2; % 设置边框的线宽
        
        % % set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.32]);
        % % set (gca,'position',[0.15,0.11,0.8,0.8] );

        figname = string(['/Users/liq/research/papers/paper15_Vesicles/01_CAC_LagrangeMultiplier/paper_Vesicle_00_Lagrange_CAC/figure_Vesicle_00_Lagrange_CAC/', pde.name, '_energy.png']);
        print(figname,'-dpng', '-r300')
    end
end

%% Mass Difference
for i_2 = 1:length(scheme2_array)
    for i_1 = 1:length(scheme1_array)

        figure(index_fig);
        index_fig = index_fig + 1;

        scheme1 = scheme1_array{i_1};
        scheme2 = scheme2_array{i_2};

        scheme = [scheme0,scheme1,scheme2];
        para.name = [pdename,[scheme1,scheme2]];

        pde = ex13_Vesicles_data_eight(para);
        hold on;
        for k = 1:maxIt
            figname_mass = [pde.name,'_S1_',num2str(pde.S1),'_dt_',num2str(dt_array(k)),'_mass.txt'];

            mass=load(figname_mass);
            tmp = 1;
            plot(mass(tmp:1:end,1),(mass(tmp:1:end,2)-mass(1,2))./mass(1,2),'-','LineWidth',2.5);
        end
        h = legend(legend_str);
        xlabel('Time','Fontsize',figure_FontSize,'interpreter','latex');
        ylabel('Ratio of Mass Difference','Fontsize',figure_FontSize-6,'interpreter','latex');
        set(gca,'FontSize',figure_FontSize-6);
        % set(gca,'linewidth',1.5);
        % % xlim([0 18])
        %         ylim([2 16])
        set(h,'box','off','interpreter','latex','FontSize',figure_FontSize-8);
        set(gca,'XTick',0:0.5:2)
        % set(gca,'YTick',-25:5:10)
        % % annotation('arrow',[0.435,0.35],[0.44,0.39],'LineWidth',2)
        box on;
        grid on;

        % 获取当前坐标轴
        ax = gca;

        % 设置坐标轴的边框粗细
        ax.Box = 'on'; % 确保边框打开
        ax.LineWidth = 2; % 设置边框的线宽

        % % set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.32]);
        % % set (gca,'position',[0.15,0.11,0.8,0.8] );

        figname = string(['/Users/liq/research/papers/paper15_Vesicles/01_CAC_LagrangeMultiplier/paper_Vesicle_00_Lagrange_CAC/figure_Vesicle_00_Lagrange_CAC/', pde.name, '_mass.png']);
        print(figname,'-dpng', '-r300')
    end
end

%% surface Area
for i_2 = 1:length(scheme2_array)
    for i_1 = 1:length(scheme1_array)

        figure(index_fig);
        index_fig = index_fig + 1;

        scheme1 = scheme1_array{i_1};
        scheme2 = scheme2_array{i_2};

        scheme = [scheme0,scheme1,scheme2];
        para.name = [pdename,[scheme1,scheme2]];

        pde = ex13_Vesicles_data_eight(para);
        hold on;
        for k = 1:maxIt
            figname_mass = [pde.name,'_S1_',num2str(pde.S1),'_dt_',num2str(dt_array(k)),'_mass.txt'];

            mass=load(figname_mass);
            tmp = 1;
            plot(mass(tmp:1:end,1),((mass(tmp:1:end,3)-mass(1,3))./mass(1,3)),'-','LineWidth',2.5);
        end
        h = legend(legend_str);
        xlabel('Time','Fontsize',figure_FontSize,'interpreter','latex');
        ylabel('Ratio of Surface Area Difference','Fontsize',figure_FontSize-6,'interpreter','latex');
        set(gca,'FontSize',figure_FontSize-6);
        % set(gca,'linewidth',1.5);
        % % xlim([0 18])
        if 1 == strcmp(scheme1,'LM4_EX_LM_')
            ylim([-4e-5 8e-5])
            %             ylim([-4e-6 8e-6])
        else
            ylim([-4e-4 8e-4])
        end
        set(h,'box','off','interpreter','latex','FontSize',figure_FontSize-8);
        set(gca,'XTick',0:0.5:2)
        % set(gca,'YTick',-25:5:10)
        % % annotation('arrow',[0.435,0.35],[0.44,0.39],'LineWidth',2)
        box on;
        grid on

        % 获取当前坐标轴
        ax = gca;

        % 设置坐标轴的边框粗细
        ax.Box = 'on'; % 确保边框打开
        ax.LineWidth = 2; % 设置边框的线宽

        % % set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.32]);
        % % set (gca,'position',[0.15,0.11,0.8,0.8] );

        figname = string(['/Users/liq/research/papers/paper15_Vesicles/01_CAC_LagrangeMultiplier/paper_Vesicle_00_Lagrange_CAC/figure_Vesicle_00_Lagrange_CAC/', pde.name, '_surface_area.png']);
        print(figname,'-dpng', '-r300')
    end
end




