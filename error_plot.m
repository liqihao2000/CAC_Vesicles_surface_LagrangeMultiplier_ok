% Error plot
clear; clc;
close all;
format long;
   
dirname={ 'example02_2D_Time_reference_LM0_SAV_1st_ok/'
          'example08_2D_Time_reference_LM0_SAV_BDF_ok/'
          'example02_2D_Time_reference_LM3_LM_1st_ok/'
          'example08_2D_Time_reference_LM3_LM_BDF_ok/'
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

% M_array = [0 50];
% M_array = [0 50 200 500];
M_array = [2];

n=0;
for kkk = 1:length(M_array)
    figure(kkk);
    lineType = {'-','.-','S--','*--','*-','--',':'};
    index_lineType=1;
    for kk = 1:size(dirname,1)
        para.M2 = M_array(kkk);
        name=[char(dirname(kk)),'phi_e',num2str(para.epsilon),...
              'M',num2str(para.M2),'Nx=',num2str(N),'Ny=',num2str(N)];
        A = readtable([name,'.txt']);
        dt_array  = A.dt_array;
        error_phi = A.error;

        if 1 == kk
            ax=loglog(dt_array,dt_array.^2*5000,'k:','linewidth',3);
            hold on;
            grid on;
            set(gca,'XMinorGrid','on','YMinorGrid','off','XMinorTick','on','YMinorTick','off');
        end
   
        if 1 == kk || 2 == kk
            linewidth = 6;
        else
            linewidth = 3;
        end
        loglog(dt_array,error_phi(:,1),lineType{index_lineType}, 'markersize',15,'linewidth',linewidth);
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
    legend({'$\mathcal{O}(\delta t^2)$', ...
        'Error\_$\phi$: SAV-1st', 'Error\_$\phi$: SAV-BDF', ...
        'Error\_$\phi$: LM-1st', 'Error\_$\phi$: LM-BDF'}, ...
           'Interpreter','latex','Location','southeast','Fontsize',15);
    xlabel('Time step $\delta t$','Interpreter','latex');
    ylabel('$L^2$ error','Interpreter','latex');
    set(gca,'FontSize',22);
    set(gca,'linewidth',1.1)
    
    figure_FontSize=35;
    set(get(gca,'XLabel'),'FontSize',figure_FontSize);
    set(get(gca,'YLabel'),'FontSize',figure_FontSize);
    set(findobj('FontSize',10),'FontSize',figure_FontSize); %这4句是将字体大小改为8号字，在小图里很清晰
    
    figname1 = ['/Users/liq/work/10_UndergraduateThesis/2024/2024_chenhuiyi/CHD_Bachelor_ChenHuiyi/figs/error_',num2str(n),'.png'];
%     print(figname1,'-dpng', '-r300')
    
end

