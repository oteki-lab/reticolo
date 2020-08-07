close all

Mx=25;
npoints=400;
h2=0.010;
diameter_x=0.100;
kappa=0.01;

periodicity_x=0.400;

h1=0.100-h2;

load('peaks_AbsRCWA_wav300_920kappa0.01_diam100_per400_h1_90_h2_10_Fourier25_nbpoints400.mat');

for i=1
    load(['FieldRCWA_','wav',num2str(locs(i)*1000),'kappa',num2str(kappa),'_diam',int2str(diameter_x*1000),'_per',int2str(periodicity_x*1000),'_h1_',int2str(h1*1000),'_h2_',int2str(h2*1000),'_Fourier',int2str(Mx),'_nbpoints1.mat']);
    
    figure
    pcolor(xx,zz,abs(Ex).^2);shading interp;colorbar;colormap(jet);hold on;
    contour(xx,zz,real(indice),'w','linewidth',1.5)
    set(gca,'clim',[0 25]);
    print(['Resonance ' num2str(26-i) 'Ex'], '-dpng')
    
    figure
    pcolor(xx,zz,abs(Ey).^2);shading interp;colorbar;colormap(jet);hold on;
    contour(xx,zz,real(indice),'w','linewidth',1.5)
    set(gca,'clim',[0 25]);
    print(['Resonance ' num2str(26-i) 'Ey'], '-dpng')
    
    figure
    pcolor(xx,zz,abs(Ez).^2);shading interp;colorbar;colormap(jet);hold on;
    contour(xx,zz,real(indice),'w','linewidth',1.5);
    set(gca,'clim',[0 10]);
    print(['Resonance ' num2str(26-i) 'Ez'], '-dpng')
end