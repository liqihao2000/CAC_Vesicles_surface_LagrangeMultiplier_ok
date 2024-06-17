clear;
clc;

dirname = 'ex16_BDF2_Vesicles_data_LM1';

datadir = [dirname,'/data'];
figdir  = [dirname,'/',dirname];

X = load([datadir, '/X.txt']);
Y = load([datadir, '/Y.txt']);

figure(5)
for t= 0:0.004:0.2
    
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


