% Error plot
clear; clc;
close all;
format long;

dirname={ 'example02_3D_Time_reference_LM0_SAV_1st_ok/'
    'example08_3D_Time_reference_LM0_SAV_BDF_ok/'
    'example02_3D_Time_reference_LM3_LM_1st_ok/'
    'example08_3D_Time_reference_LM3_LM_BDF_ok/'
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
    
    lineType = {'-','.-','S--','*--','*-','--',':'};
    index_lineType=1;
    for kk = 1:size(dirname,1)
        para.M = 2;
        name=[char(dirname(kk)),'phi_e',num2str(para.epsilon),...
            'M',num2str(para.M),'S1=',num2str(para.S1),'Nx=',num2str(N),'Ny=',num2str(N),'Nz=',num2str(N)];
        A = readtable([name,'.txt']);
        dt_array  = A.dt_array;
        error_phi = A.error;

        if 1 == kk
            ax=loglog(dt_array,dt_array*20,'k-.','linewidth',3);
            hold on;
            grid on;
            ax=loglog(dt_array,dt_array.^2*16000,'k:','linewidth',3);
            set(gca,'XMinorGrid','on','YMinorGrid','off','XMinorTick','on','YMinorTick','off');
        end

        if 1 == kk || 2 == kk
            linewidth = 6;
        else
            linewidth = 3;
        end
        loglog(dt_array,error_phi(:,1),lineType{index_lineType}, 'markersize',15,'linewidth',linewidth);

        ylim([10e-10 10e-3])

        index_lineType = index_lineType + 1;


        %     loglog(dt_array,error_phi(:,1),'S-', 'markersize',12,'linewidth',3);
        %     loglog(dt_array,error_uv(:,1),'<-', 'markersize',12,'linewidth',3);
        %     loglog(dt_array,error_u(:,1),'.-', 'markersize',30,'linewidth',3);
        %     loglog(dt_array,error_p(:,1),'*-', 'markersize',10,'linewidth',3);
        %     loglog(dt_array,erroru(:,3),'--', 'markersize',30,'linewidth',3.5);
        %     loglog(dt_array,errorv(:,3),':', 'markersize',10,'linewidth',3)


        %     xlim([10e-4/4 10e-3])
        %     if n==1
        %         ylim([10e-11 10e-5/1])
        %     end
        %     if n==2
        %         ylim([10e-11 10e0*0.2])
        %     end
    end
    legend({'$\mathcal{O}(\delta t)$', '$\mathcal{O}(\delta t^2)$',...
        'Error\_$\phi$: SAV-1st', 'Error\_$\phi$: SAV-BDF', ...
        'Error\_$\phi$: LM-1st', 'Error\_$\phi$: LM-BDF'}, ...
        'Interpreter','latex','Location','southeast','Fontsize',17);
    xlabel('Time step $\delta t$','Interpreter','latex');
    ylabel('$L^2$ error','Interpreter','latex');
    set(gca,'FontSize',22);
    set(gca,'linewidth',1.1)

    figure_FontSize=35;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize);
    set(get(gca,'YLabel'),'FontSize',figure_FontSize);
    set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰

    figname1 = ['/Users/liq/work/05_AcademicActivity/2024_05_兰州大学/Vesicle_slide/figure_vesicle_01_Lagrange/error_',num2str(n),'.png'];
    print(figname1,'-dpng', '-r300')

    n = n+1;

end

