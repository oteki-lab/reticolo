
%% Exemples ? mettre en ligne 273 du code [dans la boucle for az=1:Nb_couches, ligne avant a{az}=retcouche(init,u{az},op_retcouche);]
%% Pour d�finir diff�rentes structures

%%%%%%%%%%%% Principe de d�finition des structures ? l'int�rieur d'une p�riode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Si la couche est structur�e (Nm non-nul), on d�finit une "inclusion" d'indice Nm dans une "matrice" d'indice N
%% L'inclusion est un rectangle (Ntre=1) ou une ellipse (Ntre=5) de taille diameter_x en x et diameter_y en y
%% Les 2 premiers �l�ments du vecteur donnent la position en (x,y) du centre de cette inclusion

%% Pour faire des formes diff�rentes, on peut d�finir plusieurs inclusions, voir les exemples d'une croix et d'un L simplement form�s de 2 rectangles
%% Il est donc �ventuellement possible de faire une matrice de "pixels" enconsid�rant chaque pixel comme une inclusion...
%% Mais il est beaucoup plus facile de d�finir des formes simples par quelques d'inclusions (2 pour un L ou un +, 3 pour un H ou un U...)
%% ATTENTION : les inclusions peuvent se chevaucher, mais elles se "recouvrent" les unes les autres. 
%% C'est la derni�re qui d�finira la valeur de l'indice.


%% Quelques exemples

periodicity_x=0.35;                   % p�riode en x
periodicity_y=periodicity_x;          % p�riode en y
diameter_x=0.14;                       
diameter_y=diameter_x;                  

%% Rectangle centr? en (0,0) pour Ntre=1 ou ellipse pour Ntre=5
%% Carr? et cercle si diameter_x=diameter_y
if Nm(az)==0;u{az}=retu(period,{N(az),k0});else u{az}=retu(period,{N(az),[0,0,diameter_x,diameter_y,Nm(az),Ntre],k0});end;    

%% Cross centered in (0,0)
if Nm(az)==0;u{az}=retu(period,{N(az),k0});else u{az}=retu(period,{N(az),[0,0,diameter_x,diameter_y/3,Nm(az),1],[0,0,diameter_x/3,diameter_y,Nm(az),1],k0});end;

%% Non-centered L-shape
if Nm(az)==0;u{az}=retu(period,{N(az),k0});else u{az}=retu(period,{N(az),[-diameter_x/6,0,diameter_x/3,diameter_y,Nm(az),1],[diameter_x/2,-diameter_y/4,diameter_x,diameter_y/2,Nm(az),1],k0});end;

%% i centr? (pour se faire plaisir...)
if Nm(az)==0;u{az}=retu(period,{N(az),k0});else u{az}=retu(period,{N(az),[0,0,diameter_x/4,diameter_y,Nm(az),1],[0,diameter_y/2+diameter_y/4+diameter_y/8,diameter_x/4,diameter_y/4,Nm(az),5],k0});end;
