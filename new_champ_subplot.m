function [] = new_champ_subplot(~)

%%%% Exemple de figure pour tracer le champ pour une calcul fait avec trace_champ=1, x0=[], y0=0, z0=[]
%%%% Hy en fonction de x et z
%%%% Les contours de l'objet sont en blancs (indice contient la carte des indices) 
% figure
% pcolor(xx,zz,abs(Hy).^2);shading flat;colorbar;colormap(hot);hold on;
% contour(xx,zz,indice,'w','linewidth',1.5)


close all


[FileName,PathName] = uigetfile('*.mat','sélectionnez le fichier .mat');
 load(fullfile(PathName,FileName));
 
%Pour le tracé des cartes de champ 
 
a1=0;
a2=15;         % Definit la valeur de l'intensité maximale sur la figure pour eviter les effets de pointe
%b1=0.04;    %0.04
%b2=0.55;    %0.55 pour 200nm,0.4 pour 100nm?       % Agit sur l'extension verticale de l'image dans la profondeur de la cellule (pour voir tout le motif par exemple)
dossier='/local/home/nicolas.vandamme/Documents/MATLAB/RETICOLO2/Champs/hotcarrier';
FileName2=regexprep(FileName,'[_]','-');


% Retrouver la coupe

dzeta=0;
if isempty(x0)
    if isempty(y0)
        coupe=['z=',int2str(z0*1000)];
        dzeta=3;
    else
        coupe=['y=',int2str(y0*1000)];
        dzeta=2;
    end
else
    coupe=['x=',int2str(x0*1000)];
    dzeta=1;
end

if dzeta==0
    error('Pas de coupe de champ definie...')
end

% retrouver la polarisationh

if pol==2
        pola='TM';pola2='(H along y)';
    elseif pol==0
        pola='TE';pola2='(E along y)';
    else
        pola='undef';
end


%%%%%
switch dzeta
    
    case 2
%si coupe en y, les variables sont xx, zz

E2=figure;
titre1=['E^2 at ',int2str(wavelength*1000),'nm - coupe ',coupe,' in ',pola,' polarization ',pola2];
pcolor(xx,zz,abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);
shading interp;
colormap(hot);
title(titre1,'FontSize',15);
hold on;
contour(xx,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('z')

caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

savefigchamp2(E2,titre1,dossier)
oldFolder=cd(dossier);
    saveas(E2,titre1,'fig') 
cd(oldFolder);

%%%%%

plotall=figure('units','normalized','outerposition',[0 0 1 1]);
supertitle=['Field intensity cross section at ',int2str(wavelength*1000),'nm - coupe ',coupe,' in ',pola,' polarization ',pola2];


subplot(2,3,1)
titrea=['Ex^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,zz,abs(Ex).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('z')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,2)
titrea=['Ey^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,zz,abs(Ey).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('z')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,3)
titrea=['Ez^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,zz,abs(Ez).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('z')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,4)
titrea=['Hx^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,zz,abs(Hx).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('z')
%caxis([a1 a2]); No c axis for H to be able to compare data between curves ? 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,5)
titrea=['Hy^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,zz,abs(Hy).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('z')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,6)
titrea=['Hz^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,zz,abs(Hz).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('z')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

suplabel(['Maps obtained with calculation:',FileName2]);
[~,h3]=suplabel(supertitle,'t');
set(h3,'FontSize',30) 

savefigchamp2(plotall,supertitle,dossier)
oldFolder=cd(dossier);
    saveas(plotall,supertitle,'fig') 
cd(oldFolder);

    case 1
%si coupe en x, les variables sont yy, zz

E2=figure;
titre1=['E^2 at ',int2str(wavelength*1000),'nm - coupe ',coupe,' in ',pola,' polarization ',pola2];
pcolor(yy,zz,abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);
shading interp;
colormap(hot);
title(titre1,'FontSize',15);
hold on;
contour(yy,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('y')
ylabel('z')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

savefigchamp2(E2,titre1,dossier)
oldFolder=cd(dossier);
    saveas(E2,titre1,'fig') 
cd(oldFolder);

%%%%%

plotall=figure('units','normalized','outerposition',[0 0 1 1]);
supertitle=['Field intensity cross section at ',int2str(wavelength*1000),'nm - coupe ',coupe,' in ',pola,' polarization ',pola2];


subplot(2,3,1)
titrea=['Ex^2 at ',int2str(wavelength*1000),'nm'];
pcolor(yy,zz,abs(Ex).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(yy,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('y')
ylabel('z')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,2)
titrea=['Ey^2 at ',int2str(wavelength*1000),'nm'];
pcolor(yy,zz,abs(Ey).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(yy,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('y')
ylabel('z')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,3)
titrea=['Ez^2 at ',int2str(wavelength*1000),'nm'];
pcolor(yy,zz,abs(Ez).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(yy,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('y')
ylabel('z')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,4)
titrea=['Hx^2 at ',int2str(wavelength*1000),'nm'];
pcolor(yy,zz,abs(Hx).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(yy,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('y')
ylabel('z')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,5)
titrea=['Hy^2 at ',int2str(wavelength*1000),'nm'];
pcolor(yy,zz,abs(Hy).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(yy,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('y')
ylabel('z')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,6)
titrea=['Hz^2 at ',int2str(wavelength*1000),'nm'];
pcolor(yy,zz,abs(Hz).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(yy,zz,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [0:0.1:max(zz)]);
%set(gca, 'YTicklabel', [])
xlabel('y')
ylabel('z')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

suplabel(['Maps obtained with calculation:',FileName2]);
[~,h3]=suplabel(supertitle,'t');
set(h3,'FontSize',30) 

savefigchamp2(plotall,supertitle,dossier)
oldFolder=cd(dossier);
    saveas(plotall,supertitle,'fig') 
cd(oldFolder);

    case 3
        
       %si coupe en y, les variables sont xx, zz

E2=figure;
titre1=['E^2 at ',int2str(wavelength*1000),'nm - coupe ',coupe,' in ',pola,' polarization ',pola2];
pcolor(xx,yy,abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);
shading interp;
colormap(hot);
title(titre1,'FontSize',15);
hold on;
contour(xx,yy,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('y')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

savefigchamp2(E2,titre1,dossier)
oldFolder=cd(dossier);
    saveas(E2,titre1,'fig') 
cd(oldFolder);

%%%%%

plotall=figure('units','normalized','outerposition',[0 0 1 1]);
supertitle=['Field intensity cross section at ',int2str(wavelength*1000),'nm - coupe ',coupe,' in ',pola,' polarization ',pola2];


subplot(2,3,1)
titrea=['Ex^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,yy,abs(Ex).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,yy,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('y')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,2)
titrea=['Ey^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,yy,abs(Ey).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,yy,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('y')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,3)
titrea=['Ez^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,yy,abs(Ez).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,yy,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('y')
caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,4)
titrea=['Hx^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,yy,abs(Hx).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,yy,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('y')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,5)
titrea=['Hy^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,yy,abs(Hy).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,yy,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('y')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

%%%%%
subplot(2,3,6)
titrea=['Hz^2 at ',int2str(wavelength*1000),'nm'];
pcolor(xx,yy,abs(Hz).^2);
shading interp;
colormap(hot);
title(titrea,'FontSize',15);
hold on;
contour(xx,yy,indice,'w','linewidth',1.5)

set(gca, 'XTick', [-(periodicity_x/2):0.1:(periodicity_x/2)]);
%set(gca, 'XTicklabel', [])
set(gca, 'YTick', [-(periodicity_y/2):0.1:(periodicity_y/2)]);
%set(gca, 'YTicklabel', [])
xlabel('x')
ylabel('y')
%caxis([a1 a2]); 
set(gcf,'Renderer','Zbuffer');
colorbar;

suplabel(['Maps obtained with calculation:',FileName2]);
[~,h3]=suplabel(supertitle,'t');
set(h3,'FontSize',30) 

savefigchamp2(plotall,supertitle,dossier)
oldFolder=cd(dossier);
    saveas(plotall,supertitle,'fig') 
cd(oldFolder); 

end

close all


end

% syntax : savefig(figname,filename)
% sauve dans un dossier (nomme Images cree dans le dossier courant) la figure intitulee 'figname' sous le nom
% filename en eps integrable dans Latex. La resolution de l'image est
% soi-disant celle de l'ecran.(-r0)
function savefigchamp2(figname,filename,folder)
Home = pwd;
% if nargin<3
%     Dossier = dir('Images');
%         if size(Dossier,1) == 0;
%             mkdir ('Images')
%         end
%     cd Images
% else

    if (isdir(folder)==1)
        cd(folder)
    else
        mkdir(folder)
        cd(folder)
    end
% end
% Mise en current figure
set(0,'CurrentFigure',figname);
% Mise en forme pour impression en pdf
set(gcf,'PaperType','A5')
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize',[23 18])
set(gcf, 'PaperPosition', [-4 -3 30 23]);
%print ([ '-f' num2str(figname)], '-depsc2','-r100',filename);
print ([ '-f' num2str(figname)], '-dpdf','-r200','-loose',filename);
print ([ '-f' num2str(figname)], '-dpng','-r200','-loose',filename);


cd (Home)
end