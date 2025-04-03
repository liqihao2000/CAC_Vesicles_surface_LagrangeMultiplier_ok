% Error plot
clear; clc;
close all;
format long;

dirname={ 
    'example02_3D_Time_reference_MSAV_1st_ok/'
%     'example02_3D_Time_reference_LM0_SAV_1st_ok/'
%     'example02_3D_Time_reference_LM1_SAV_1st_ok/'
%     'example02_3D_Time_reference_LM3_LM_1st_ok/'
    'example08_3D_Time_reference_MSAV_BDF_ok/'

    'example02_3D_Time_reference_LM4_EX_LM_1st_ok/'

%     'example08_3D_Time_reference_MSAV_BDF_ok/'
%     'example08_3D_Time_reference_LM0_SAV_BDF_ok/'
%     'example08_3D_Time_reference_LM1_SAV_BDF_ok/'
%     'example08_3D_Time_reference_LM3_LM_BDF_ok/'
    'example08_3D_Time_reference_LM4_EX_LM_BDF_ok/'
    };

% Space: Domain and N
domain.left   = 0;
domain.right  = 2*pi;
domain.bottom = 0;
domain.top    = 2*pi;

Lx = domain.right - domain.left;
Ly = domain.top   - domain.bottom;

% Parameters
para.epsilon = 6*pi/128;
para.M = 2;
N=32;

M = 2;
newcolors = [
        0   0.447000000000000   0.741000000000000
        0.850000000000000   0.325000000000000   0.098000000000000
        0.929000000000000   0.694000000000000   0.125000000000000
        0.301000000000000   0.745000000000000   0.933000000000000
        0.266000000000000   0.574000000000000   0.188000000000000
        0.494000000000000   0.184000000000000   0.556000000000000
        0.635000000000000   0.078000000000000   0.184000000000000
        ];

n=3;
for kkk = 1:2
    if 1 == kkk
        para.S1 = 0;
        para.S2 = 0;
        para.S3 = 0;
    elseif 2 == kkk
        para.S1 = 4;
        para.S2 = 4;
        para.S3 = 1;
    end

    figure(kkk+2);
    
    colororder(newcolors)
    
    lineType = {'*-','.-','S--','^--','*-','--',':'};
    index_lineType=1;
    for kk = 1:size(dirname,1)
        para.M = 2;
        name=[char(dirname(kk)),'phi_e',num2str(para.epsilon),...
            'M',num2str(para.M),'S1=',num2str(para.S1),'Nx=',num2str(N),'Ny=',num2str(N),'Nz=',num2str(N)];
        A = readtable([name,'.txt']);
        dt_array  = A.dt_array;
        error_phi = A.error;

        if 1 == kk
            loglog(dt_array,dt_array*20,'k-.','linewidth',3);
            hold on;
            grid on;
            ax=loglog(dt_array,dt_array.^2*16000,'k:','linewidth',3,'MarkerSize',30);
            set(gca,'XMinorGrid','on','YMinorGrid','off','XMinorTick','on','YMinorTick','off');
        end

        if  2 == kk
            linewidth = 6;
            markersize = 60;
        else
            linewidth = 4;
            markersize = 20;
        end
        loglog(dt_array,error_phi(:,1),lineType{index_lineType}, 'markersize',markersize,'linewidth',linewidth);

        ylim([10e-10 10e-1])

        index_lineType = index_lineType + 1;
    end
%     legend({'$\mathcal{O}(\delta t)$', '$\mathcal{O}(\delta t^2)$',...
%         'Error\_$\phi$: MSAV-1st', 'Error\_$\phi$: MSAV-BDF', ...
%         'Error\_$\phi$: CSAV-1st', 'Error\_$\phi$: CSAV-BDF'}, ...
%         'Interpreter','latex','Location','southeast','Fontsize',17);
    legend({'$\mathcal{O}(\delta t)$', '$\mathcal{O}(\delta t^2)$',...
        'MSAV-1st', 'MSAV-BDF', 'CSAV-1st', 'CSAV-BDF'}, ...
        'Interpreter','latex','Location','southeast','Fontsize',17);
    xlabel('Time step $\delta t$','Interpreter','latex');
    ylabel('$L^2$ error','Interpreter','latex');
    set(gca,'FontSize',22);
    set(gca,'linewidth',2);
%     axis square;
%     daspect([1.1 10000 20000])
%     d = daspect

%       pbaspect([1 2 1]); % 设置纵横比为 1:2
%     daspect([2000, 1, 1]);
%     set(gca, 'DataAspectRatio', [1, 1, 1]);
%     set(gcf, 'PaperPosition', [0 0 80 6]);
%     set(gca, 'Position', [0.22, 0.18, 0.6, 0.84]); % 修改窗口大小和位置
%     ax = gca; % 获取当前坐标轴句柄
%     ax.DataAspectRatio = [newWidthRatio, 1, 1]; 

    set(gcf, 'Position', [100, 100, 560, 470]); % 调整图形窗口大小
    
    figure_FontSize=30;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize);
    set(get(gca,'YLabel'),'FontSize',figure_FontSize);
    set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰

    figname1 = ['/Users/liq/research/papers/paper15_Vesicles/01_CAC_LagrangeMultiplier/paper_Vesicle_00_Lagrange_CAC/figure_Vesicle_00_Lagrange_CAC/error_',num2str(n),'.png'];
%     print(figname1,'-dpng', '-r300')

    n = n+1;

end


