clear;
clc;

dirname = 'ex15_1st_Vesicles_data_LM_1st_old';
% dirname = 'ex15_1st_Vesicles_data_LM_bdf_old';

dirname = 'ex15_LM0_SAV_1st_Vesicles';
dirname = 'ex15_LM0_SAV_BDF_Vesicles';
dirname = 'ex15_LM1_SAV_1st_Vesicles';
dirname = 'ex15_LM1_SAV_BDF_Vesicles';
dirname = 'ex15_LM3_LM_1st_Vesicles';
dirname = 'ex15_LM3_LM_BDF_Vesicles';


datadir = [dirname,'/data'];
figdir  = [dirname,'/',dirname];

X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);

figure(5)
for t= 0:0.4:20
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


