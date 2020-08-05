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