function []=parfor_cellule_GaAs_ISE(FF,hauteur_struc)
%%%%                Calcul de la diffraction par un r�seau 2D compos� de 1 � 12 couches
%%%%                                         (17/06/2014)
%%%%                                    
%%%% Chaque couche peut �tre structur�e ou non
%%%% Toutes les couches ont la m�me p�riode (rectangulaire) et les m�mes motifs (rectangulaires)
%%%% Le r�seau est suppos� sub-longueur d'onde (p�riode < longueur d'onde), seul l'ordre (0,0) est calcul� dans l'air
%%%% L'empilement vertical se termine par un substrat m�tallique (en tout cas absorbant, on ne calcule pas la transmission)
%%%% 
%%%% Grandeurs calcul�es en fonction de la longueur d'onde : 
%%%% - Coefficient de r�flexion dans l'ordre 0
%%%% - Absorption dans les diff�rentes couches 
%%%% 
%%%% On peut �galement calculer pour une longueur d'onde
%%%% - Le champ �lectromagn�tique en volume dans une couche
%%%% - Une coupe du champ �lectromagn�tique (au choix parmi les 3 possibilit�s x=x0, y=y0 ou z=z0)
%%%%
%%%% On fait le calcul soit avec RCWA standard (op_granet=0) soit avec RCWA modifi�e (op_granet=1) pour am�liorer la convergence
%%%% Cette option peut �tre utile pour les r�seaux m�talliques. A priori pas n�cessaire pour les r�seaux di�lectriques
%%%%
%%%% Unit� de longueur : �m  (la longueur d'onde et toutes les longueurs sont en �m)
%%%%
%
% Empilement (vue de c�t� de la structure)
% -----------------------------------------
% Les couches sont num�rot�es i : indices ni et nim ; hauteur hi.
% nim est l'indice dans les plots (cf. vue de dessus) et ni est l'indice entre les plots
% ATTENTION : Si la couche n'est pas structur�e, mettre nim=0
% Le nombre de couches est NB_couches, elles sont num�rot�es de 1 � 12 en partant du haut
% S'IL Y A UNE SEULE COUCHE, C'EST LA COUCHE 1. S'IL Y A DEUX COUCHES, COUCHES 1 et 2, etc.
% ATTENTION : Pour les couches qui n'existent pas, mettre hi=0
% 
% L'indice du milieu incident est nh (milieu homog�ne transparent semi-infini)
% Le substrat est un m�tal (ou un milieu absorbant) d'indice nsub
%                                                                                                                                  
%          milieu incident (nh)                                                                                                   
%      ^ ***************************************************************************                                                                                                
%   h1 |    n1      *  n1m    *    n1    *    n1m  *    n1    *   n1m   *     n1                                    
%      x ***************************************************************************
%   h2 |    n2      *  n2m    *    n2    *    n2m  *    n2    *   n2m   *     n2                              
%      x ***************************************************************************                                                              
%                                  .  
%                                  .
%                                  .
%      x ***************************************************************************                                                                                            
%   h6 |    n6      *  n6m    *    n6    *    n6m  *    n6    *   n6m   *     n6                                                                                                                                                                                                                                       
%      x ***************************************************************************                                                                                                                                                                                                            
%           substrat m�tallique (nsub)                                                                                                       
% z ^                                                                                   
%   |
%   --->  x
%
% Vue de dessus de la couche i
% -----------------------------
%                              periodicity_x
%                             <------------->                                                                             
%         *******       *******       *******       *******       ******* ^
%         *******       *******       *******       *******       ******* |                                                                  
%         *******       *******       *******       *******       ******* | periodicity_y                                                             
%                                                                         |                                                                                                                                         
%                                                                         |                                
%         *******       *******       *******       *******       ******* v                                                                   
%         *******       **nim**       *******       *******       *******                                                                    
%         *******       *******       *******       *******       *******
%                                ni                                                                                                                                                                                  
%                                                                                                         
%         *******       *******       *******       *******       ******* ^                                                                   
% y ^     *******       *******       *******       *******       ******* | diameter_y                                                                  
%   |     *******       *******       *******       *******       ******* v                                                                  
%   |                                 <----->
%   --->  x                           diameter_x
%                                                                                                                                                                                         
%
% Sorties du code (vecteurs de la m�me taille que wavelength ou matrices de taills Nb_couches x length(wavelength), sauf les champs)
% ------------------------------------------------------------------------------------------------------------------------------------
% On calcule les coefficients de r�flexion dans l'ordre (0,0) uniquement (r�seau suppos� sub-longueur d'onde, ordres sup�rieurs �vanescents)
%
% ref_TE_TE : coeff de r�flexion en amplitude pour une onde incidente TE et une onde r�fl�chie TE
% ref_TE_TM : coeff de r�flexion en amplitude pour une onde incidente TE et une onde r�fl�chie TM
% ref_TM_TM : coeff de r�flexion en amplitude pour une onde incidente TM et une onde r�fl�chie TM
% ref_TM_TE : coeff de r�flexion en amplitude pour une onde incidente TM et une onde r�fl�chie TE
% R0_TE_TE : coeff de r�flexion en intensit� pour TE --> TE 
% R0_TE_TM : coeff de r�flexion en intensit� pour TE --> TM 
% R0_TM_TM : coeff de r�flexion en intensit� pour TM --> TM
% R0_TM_TE : coeff de r�flexion en intensit� pour TM --> TE
%
% Note : Il n'y a que quand theta(1)~=0 ET theta(2)~=0 que les 4 coefficients de r�flexion sont non nuls. 
% Sinon, d�s que l'un des angles theta est �gal � z�ro, le programme utilise la symm�trie du probl�me et un seul coefficient est non nul.
% La polarisation de l'onde incidente est dans ce cas fix�e avec la variable pol (cf param�tres num�riques)
%
% Si cal_abs=1
% --------------
% Abs_plots(i,:) : absorption dans les plots de la couche i en fonction de la longueur d'onde
% Abs(i,:) : absorption dans toute la couche i (dans les plots + entre les plots) en fonction de la longueur d'onde
% A_sub : absorption dans le substrat m�tallique 
% A_tot : absorption totale dans tout l'empilement
% test : 1-R0-A_tot. Doit �tre petit
%
% Si cal_champ=1  (dans ce cas, une seule longueur d'onde)
% ----------------------------------------------------------
% E_semicon, H_semicon : composantes des champs E et H dans la couche dont le num�ro a �t� sp�cifi� par la variable N_semicon
%                        le champ est calcul� aux points (x_semicon,y_semicon,z_semicon)
%                        vecteurs de dimensions (Nb_pts_z_semicon,165,150,3), points en z, x, y, et composantes x, y et z 
%                        E_semicon(:,:,:,1) est Ex en tous points
% W_semicon : permet de calculer l'int�grale du champ E_semicon (m�thode de Gauss)
%             int_{maille}|E|^2 = sum(sum(sum(W_semicon.*sum(abs(E_semicon).^2,4))))
% Einc, Hinc : composantes du champ E et du champ H de l'onde plane incidente utilis�e pour calculer le champ E_semicon et H_semicon
%
% Si trace_champ=1  (dans ce cas, une seule longueur d'onde)
% ------------------------------------------------------------
% Ex,Ey,Ez : coupe du champ E aux points zz,xx,yy 
% Hx,Hy,Hz : coupe du champ H aux points zz,xx,yy
% On choisit une coupe parmi x=x0, y=y0 ou z=z0

%clear;
%retio;
%[prv,vmax]=retio([],inf*i);          % PAS �criture sur fichiers (inutile sur des versions r�centes de Matlab)

%%%%%% Nombre de couches dans l'empilement
Nb_couches=6;                        % Nombre de couches dans l'empilement, chiffre compris entre 1 et 12
                                     % Les couches sont num�rot�es de haut en bas

%%%%%% Longueur d'onde et angle d'incidence
wavelength=[0.3:0.01:0.5,0.505:0.005:0.68,0.682:0.002:0.9,0.905:0.005:0.95];
theta=[0,0];                          % angle d'incidence en degr�s

%%%%%% Param�tres g�om�triques
periodicity_x=0.7;                   % p�riode en x
periodicity_y=periodicity_x;          % p�riode en y
diameter_x=FF*periode;                      % taille des plots en x 
diameter_y=diameter_x;                % taille des plots en y 
h1=0.06;                                 % �paisseur de la couche 1 (0 si pas de couche)
h2=0.025;                                 % �paisseur de la couche 2 (0 si pas de couche)
h3=0.205;                                 % �paisseur de la couche 3 (0 si pas de couche)
h4=0.100;                                 % �paisseur de la couche 4 (0 si pas de couche)
h5=hauteur_struc;                                 % �paisseur de la couche 5 (0 si pas de couche)
h6=0.1; 
h7=0;                                 % �paisseur de la couche 7 (0 si pas de couche)
h8=0;                                 % �paisseur de la couche 8 (0 si pas de couche)
h9=0;                                 % �paisseur de la couche 9 (0 si pas de couche)
h10=0;                                % �paisseur de la couche 10 (0 si pas de couche)
h11=0;                                % �paisseur de la couche 11 (0 si pas de couche)
h12=0;                                % �paisseur de la couche 12 (0 si pas de couche)

%%%%%% Indices de r�fracion (de haut en bas), �ventuellement fonction de la longueur d'onde
nh=1;                                    % indice du milieu incident
n1=retindice_compile(wavelength,23.2,'linear');     % indice entre les plots de la couche 1 
n1m=0;                                   % indice des plots de la couche 1 (0 si couche pas structur�e)
n2=retindice_compile(wavelength,4.888,'linear');       % indice entre les plots de la couche 2 
n2m=0;                                   % indice des plots de la couche 2 (0 si couche pas structur�e)
n3=retindice_compile(wavelength,4.702,'linear');     % indice entre les plots de la couche 3
n3m=0;                                   % indice des plots de la couche 3 (0 si couche pas structur�e)
n4=retindice_compile(wavelength,4.6021,'linear');     % indice entre les plots de la couche 4
n4m=0;                                   % indice des plots de la couche 4 (0 si couche pas structur�e)
n5=1.9*ones(size(wavelength));      % indice entre les plots de la couche 5
n5m=retindice_compile(wavelength,1.7,'linear');         % indice des plots de la couche 5 (0 si couche pas structur�e)
n6=retindice_compile(wavelength,1.7,'linear');    % indice entre les plots de la couche 6
n6m=0;                                   % indice des plots de la couche 6 (0 si couche pas structur�e)
n7=1;                                   % indice entre les plots de la couche 7
n7m=0;                                   % indice des plots de la couche 7 (0 si couche pas structur�e)
n8=1;                                    % indice entre les plots de la couche 8
n8m=0;                                   % indice des plots de la couche 8 (0 si couche pas structur�e)
n9=1;                                    % indice entre les plots de la couche 9
n9m=0;                                   % indice des plots de la couche 9 (0 si couche pas structur�e)
n10=1;                                   % indice entre les plots de la couche 10
n10m=0;                                  % indice des plots de la couche 10 (0 si couche pas structur�e)
n11=1;                                   % indice entre les plots de la couche 11
n11m=0;                                  % indice des plots de la couche 11 (0 si couche pas structur�e)
n12=1;                                   % indice entre les plots de la couche 12
n12m=0;                                  % indice des plots de la couche 12 (0 si couche pas structur�e)
nsub=retindice_compile(wavelength,1.7,'linear');    % indice du m�tal du substrat 

%%%%%% Param�tres num�riques
pol=2;                               % polarisation incidente, TM pol=2  TE pol=0. TE (TM) = champ E (champ H) Transverse au plan d'incidence
                                     % En incidence normale, TM veut dire H//y et TE est E//y
Mx=25;                               % nombre de termes de Fourier en x
My=Mx;                                % nombre de termes de Fourier en y
op_granet=1;                         % si 1, RCWA modifi�e pour am�liorer la convergence (transform�e de coordonn�es r�elles sur les discontinuit�s) 
Nb_pts_z=20;                         % nombre de points en z pour calculer l'absorption (dans les couches o� on calcule le champ en volume)
cal_abs=1;                           % si 1, calcul de l'absorption dans toutes les couches
                                     % si 0, seulement calcul du coeff de r�flexion
cal_champ=0;                         % si 1, calcul du champ en volume dans la couche N_semicon 
                                     % si 0, pas calcul du champ dans la couche N_semicon
                                     % ATTENTION: dans ce cas, seulement une longueur d'onde (le calcul peut �tre tr�s lourd selon Nb_pts_z_semicon)
N_semicon=3;                         % Num�ro de la couche dans laquelle on calcule le champ en volume (si cal_champ=1) 
Nb_pts_z_semicon=50;                 % nombre de points en z pour calculer le champ en volume dans la couche N_semicon (si cal_champ=1)                                        
trace_champ=0;                       % si 1, calcul d'une coupe du champ
                                     % ATTENTION: dans ce cas seulement une longueur d'onde
x0=[];                                % coupe en x=x0 si trace_champ=1 ([] si on ne fait pas de coupe en x0)
y0=0;                                % coupe en y=y0 si trace_champ=1 ([] si on ne fait pas de coupe en y0)                       
z0=[];                                % coupe en z=z0 si trace_champ=1 ([] si on ne fait pas de coupe en z0)
                                     % z=0 est situ� tout en bas de la structure, � une hauteur h_sub dans le substrat
h_air=0.05;                          % �paisseur au-dessus pour trac� du champ (trace_champ=1)
h_sub=0.05;                          % �paisseur dans le substrat pour trac� du champ (trace_champ=1)
h_2pts=20;                           % distance entre 2 pts en z pour les couches d'�paisseur sup�rieure � 1�m pour le trac� du champ (trace_champ=1)
op_objet=0;                          % si 1, trac� de la g�om�trie pour v�rifier la structure calcul�e

if op_granet==1;Bx=500;Ax=0.02/Bx;By=Bx;Ay=Ax;xdisc=[-diameter_x/2,diameter_x/2];ydisc=[-diameter_y/2,diameter_y/2];end;
% if op_granet==1;Bx=500;Ax=0.02/Bx;By=[];Ay=[];xdisc=[-diameter_x/2,diameter_x/2];ydisc=[];end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref_TE_TE=zeros(1,length(wavelength));
ref_TE_TM=ref_TE_TE;
ref_TM_TM=ref_TE_TE;
ref_TM_TE=ref_TE_TE;
R0_TE_TE=zeros(1,length(wavelength));
R0_TE_TM=R0_TE_TE;
R0_TM_TM=R0_TE_TE;
R0_TM_TE=R0_TE_TE;
A_tot=zeros(1,length(wavelength));
A_sub=zeros(1,length(wavelength));
Abs=zeros(Nb_couches,length(wavelength));
Abs_plots=zeros(Nb_couches,length(wavelength));
test=zeros(1,length(wavelength));
Einc=[];Hinc=[];E_semicon=[];H_semicon=[];x_semicon=[];y_semicon=[];z_semicon=[];W_semicon=[];
Ex=[];Ey=[];Ez=[];Hx=[];Hy=[];Hz=[];
xx=[];yy=[];zz=[];indice=[];
Ntre=1;
H=[h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12];
if cal_abs==1||cal_champ==1||trace_champ==1;op_retcouche=1;else op_retcouche=0;end;
if H(Nb_couches)<1e-5;disp('WARNING : There is a problem in the definition of the layers number !!');return;end;
if trace_champ==1&&isempty(x0)==1&&isempty(y0)==1&&isempty(z0)==1;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;
if trace_champ==1&&isempty(x0)==0&&isempty(y0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;
if trace_champ==1&&isempty(x0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;
if trace_champ==1&&isempty(y0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!');return;end;
% parfor_progress(length(wavelength))
parfor zou=1:length(wavelength)
    inc=[];
    sym=[];
    e0=[];
    o0=[];
    R0_TE_TE_vect=zeros(1);
    R0_TE_TM_vect=R0_TE_TE_vect;
    R0_TM_TM_vect=R0_TE_TE_vect;
    R0_TM_TE_vect=R0_TE_TE_vect;
    ref_TE_TE_vect=zeros(1);ref_TE_TM_vect=ref_TE_TE_vect;ref_TM_TM_vect=ref_TE_TE_vect;ref_TM_TE_vect=ref_TE_TE_vect;
    A_tot_vect=zeros(1);
    A_sub_vect=zeros(1);
    Abs_vect=zeros(Nb_couches,1);
    Abs_plots_vect=zeros(Nb_couches,1);
    test_vect=zeros(1);
    Einc=[];Hinc=[];E_semicon=[];H_semicon=[];x_semicon=[];y_semicon=[];z_semicon=[];W_semicon=[];
    Ex=[];Ey=[];Ez=[];Hx=[];Hy=[];Hz=[];
    xx=[];yy=[];zz=[];indice=[];
    Ntre=1;
    lwa=wavelength(zou);
    k0=2*pi/lwa;
    %Vstart=tic
    period=[periodicity_x,periodicity_y];

    ns=nsub(zou);
    if length(n1)==1;nn1=n1;else nn1=n1(zou);end;
    if length(n2)==1;nn2=n2;else nn2=n2(zou);end;
    if length(n3)==1;nn3=n3;else nn3=n3(zou);end;
    if length(n4)==1;nn4=n4;else nn4=n4(zou);end;
    if length(n5)==1;nn5=n5;else nn5=n5(zou);end;
    if length(n6)==1;nn6=n6;else nn6=n6(zou);end;
    if length(n7)==1;nn7=n7;else nn7=n7(zou);end;
    if length(n8)==1;nn8=n8;else nn8=n8(zou);end;
    if length(n9)==1;nn9=n9;else nn9=n9(zou);end;
    if length(n10)==1;nn10=n10;else nn10=n10(zou);end;
    if length(n11)==1;nn11=n11;else nn11=n11(zou);end;
    if length(n12)==1;nn12=n12;else nn12=n12(zou);end;
    if length(n1m)==1;nn1m=n1m;else nn1m=n1m(zou);end;
    if length(n2m)==1;nn2m=n2m;else nn2m=n2m(zou);end;
    if length(n3m)==1;nn3m=n3m;else nn3m=n3m(zou);end;
    if length(n4m)==1;nn4m=n4m;else nn4m=n4m(zou);end;
    if length(n5m)==1;nn5m=n5m;else nn5m=n5m(zou);end;
    if length(n6m)==1;nn6m=n6m;else nn6m=n6m(zou);end;
    if length(n7m)==1;nn7m=n7m;else nn7m=n7m(zou);end;
    if length(n8m)==1;nn8m=n8m;else nn8m=n8m(zou);end;
    if length(n9m)==1;nn9m=n9m;else nn9m=n9m(zou);end;
    if length(n10m)==1;nn10m=n10m;else nn10m=n10m(zou);end;
    if length(n11m)==1;nn11m=n11m;else nn11m=n11m(zou);end;
    if length(n12m)==1;nn12m=n12m;else nn12m=n12m(zou);end;
    N=[nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,nn9,nn10,nn11,nn12];
    Nm=[nn1m,nn2m,nn3m,nn4m,nn5m,nn6m,nn7m,nn8m,nn9m,nn10m,nn11m,nn12m];

    if theta(1)==0 && theta(2)==0;sym=[pol-1,pol-1,0,0];end; 
    if theta(1)~=0 && theta(2)==0;sym=[0,pol-1,0,0];end; 
    if theta(1)==0 && theta(2)~=0;sym=[1-pol,0,0,0];end;
    if theta(1)~=0 && theta(2)~=0;sym=[];end;
    kx=k0*nh*sin(theta(1)*pi/180);
    ky=k0*nh*sin(theta(2)*pi/180);
    beta=[kx,ky];
    
    uh=retu(period,{nh,k0});                                     
    ub=retu(period,{ns,k0});
    
    if op_granet==1;
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym,{[],{xdisc,Ax,Bx,ydisc,Ay,By}});ah=retcouche(init,uh,1);
    else
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym);ah=retcouche(init,uh,op_retcouche);
    end
    ab=retcouche(init,ub,op_retcouche);
    
    u=[]
    a=[]
    for az=1:Nb_couches
        if Nm(az)==0;u{az}=retu(period,{N(az),k0});else u{az}=retu(period,{N(az),[0,0,diameter_x,diameter_y,Nm(az),Ntre],k0});end;
        a{az}=retcouche(init,u{az},op_retcouche);
    end
    struct_test=[]
    if op_objet==1
        struc_test=cell(1,Nb_couches+2);
        struct_test{1}={0.1,uh};
        for az=1:Nb_couches
           struct_test{az+1}={H(az),u{az}}; 
        end
        struct_test{Nb_couches+2}={0.05,ub};
        rettestobjet(period,struct_test,[0,-1.5,1],[2,2,-period/2]);
        rettestobjet(period,struct_test,[0,-2],[2,2,-period/2]);
    end

    if op_granet==1
        init_pasgranet=retinit(period,[-Mx,Mx,-My,My],beta,sym);
        ah_pasgranet=retcouche(init_pasgranet,uh,1);
        s_passage=retpassage(init,ah_pasgranet,ah,1);
        [sh,haut,beth,ch,anglh]=retb(init_pasgranet,ah_pasgranet,1e-4);
        sh=retss(sh,s_passage,retc(ah,2*lwa));
    else
        [sh,haut,beth,ch,anglh]=retb(init,ah,1e-4);
    end
    ordre=[0,0];
    K=[kx,ky]+ordre*2*pi./period;
    inch=find(((K(1)-beth(:,1)).^2+(K(2)-beth(:,2)).^2)<1e-8);
    difh=inch;
    sh=rettronc(sh,haut(inch,1),haut(difh,1),1);

    [sb,bas,betb,cb,anglb]=retb(init,ab,-0.001);
    incb=[];difb=[];
    sb=rettronc(sb,bas(incb,1),bas(difb,1),-1);
    stemp=retc(a{1},H(1));
    for az=2:Nb_couches
        stemp=retss(stemp,retc(a{az},H(az)));
    end

    stot=retss(sh,stemp,sb);
    ef=retreseau(init,stot,betb,cb,anglb,incb,difb,beth,ch,anglh,inch,difh);
    if isempty(ef.amplitude(ef.dif.TEh,ef.inc.TEh))==0;ref_TE_TE_vect=ef.amplitude(ef.dif.TEh,ef.inc.TEh);end;
    if isempty(ef.amplitude(ef.dif.TMh,ef.inc.TEh))==0;ref_TE_TM_vect=ef.amplitude(ef.dif.TMh,ef.inc.TEh);end;
    if isempty(ef.amplitude(ef.dif.TMh,ef.inc.TMh))==0;ref_TM_TM_vect=ef.amplitude(ef.dif.TMh,ef.inc.TMh);end;
    if isempty(ef.amplitude(ef.dif.TEh,ef.inc.TMh))==0;ref_TM_TE_vect=ef.amplitude(ef.dif.TEh,ef.inc.TMh);end;
    R0_TE_TE_vect=abs(ref_TE_TE_vect)^2;R0_TE_TM_vect=abs(ref_TE_TM_vect)^2;R0_TM_TM_vect=abs(ref_TM_TM_vect)^2;R0_TM_TE_vect=abs(ref_TM_TE_vect)^2;
    R0_TE_TE(:,zou)=R0_TE_TE_vect;
    R0_TE_TM(:,zou)=R0_TE_TM_vect;
    R0_TM_TE(:,zou)=R0_TM_TE_vect;
    R0_TM_TM(:,zou)=R0_TM_TM_vect;
    ref_TE_TE(:,zou)=ref_TE_TE_vect;
    ref_TE_TM(:,zou)=ref_TE_TM_vect;
    ref_TM_TE(:,zou)=ref_TM_TE_vect;
    ref_TM_TM(:,zou)=ref_TM_TM_vect;
    if op_granet==1
        [Xdisc,Ydisc]=retgranet(init,[-diameter_x/2,diameter_x/2],[-diameter_y/2,diameter_y/2]);
        [X,wX]=retgauss(-periodicity_x/2,periodicity_x/2,15,10,Xdisc);
        [Y,wY]=retgauss(-periodicity_y/2,periodicity_y/2,15,10,Ydisc);
        [x,y]=retgranet(init,'num2phys',X,Y);
        [X,Y,Xder,Yder]=retgranet(init,x,y);
        wx=wX./Xder;wy=wY./Yder;
    else
        [x,wx]=retgauss(-periodicity_x/2,periodicity_x/2,15,10,[-diameter_x/2,diameter_x/2]);
        [y,wy]=retgauss(-periodicity_y/2,periodicity_y/2,15,10,[-diameter_y/2,diameter_y/2]);
    end
    
    if cal_abs==1||cal_champ==1||trace_champ==1
        
    sb_norm=retb(init,ah,-0.1,0,[],[]);
    tab_norm=[0,1,1];
    if size(ef.inc.teta,2)==1
        inc=1;
    elseif size(ef.inc.teta,2)==2
        inc=[0,0];
        if pol==0;inc(ef.inc.TE)=1;elseif pol==2;inc(ef.inc.TM)=1;end;
    end
    [ei,zi]=retchamp(init,{ah},sh,sb_norm,inc,{x,y},tab_norm,[],1:6,1,1,1:6);
    flux=retpoynting(ei,[0,0,-1],wx,wy,[]);
    inc=1/sqrt(flux).*inc;
    
    [einc,zinc]=retchamp(init,{ah},sh,sb_norm,inc,{0,0},tab_norm,[],1:6,1,1,1:6);
    Einc=squeeze(einc(:,:,:,1:3));
    Hinc=squeeze(einc(:,:,:,4:6));
    
    end

    if cal_abs==1
        
    nunu=find(Nm~=0);num=[];
    for zer=1:length(nunu);if nunu(zer)==1;num=[num,nunu(zer),nunu(zer)+1];else;num=[num,nunu(zer)-1,nunu(zer),nunu(zer)+1];end;end;
    num=retelimine(num);
    
    tab=[];tab2=tab;
    struct={ab};
    for az=1:Nb_couches
        tab=[tab;[H(az),az+1,0]];
        tab2=[tab2;[0,az+1,1];[H(az),az+1,0]];
        struct={struct{:},a{az}};
    end
    tab(num,3)=Nb_pts_z;
    tab2=[tab2;[0,1,1];[1,1,0];[0,1,1]];
    
    [e,z,wz,o]=retchamp(init,struct,sh,sb,inc,{x,y},tab,[],[1:6]+7.25*i,1,1,1:6);
    [e2,z2,wz2,o2]=retchamp(init,struct,sh,sb,inc,{x,y},tab2,[],[1:6]+7.25*i,1,1,1:6);
    for ii=1:3                                  
        o(:,:,:,ii+3)=o(:,:,:,ii+3)./o(:,:,:,ii);
        o(:,:,:,ii)=1;
    end
    [Wz,Wx,Wy]=ndgrid(wz,wx,wy);
    W=Wz.*Wx.*Wy;
    numx=find(x>-diameter_x/2 & x<diameter_x/2);numy=find(y>-diameter_y/2 & y<diameter_y/2);
    
    flux_poyn=retpoynting(e2,[0,0,-1],wx,wy,[]);
    
    A_sub_vect=flux_poyn(2)-flux_poyn(1);

    for az=1:Nb_couches
        Abs_vect(Nb_couches-az+1)=flux_poyn(az+2)-flux_poyn(az+1);
    end
   
    for az=1:length(num)
        Abs_plots_vect(num(az))=0.5*k0*sum(sum(sum(W((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,numx,numy).*imag(o((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,numx,numy,4)).*sum(abs(e((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,numx,numy,1:3)).^2,4))));
        Abs_vect(num(az))=0.5*k0*sum(sum(sum(W((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,:,:).*imag(o((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,:,:,4)).*sum(abs(e((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,:,:,1:3)).^2,4))));
    end
    
    Abs(:,zou)=Abs_vect;
    Abs_plots(:,zou)=Abs_plots_vect;
    A_tot_vect=flux_poyn(end)-flux_poyn(1);
    test_vect=1-R0_TE_TE_vect-R0_TE_TM_vect-R0_TM_TM_vect-R0_TM_TE_vect-A_tot_vect;
    test(:,zou)=test_vect
    A_tot(:,zou)=A_tot_vect;
    A_sub(:,zou)=A_sub_vect;
    
    end

    if cal_champ==1
    
    tab_semicon=[];
    struct={ab};
    for az=1:Nb_couches
        tab_semicon=[tab_semicon;[H(az),az+1,0]];
        struct={struct{:},a{az}};
    end
    tab_semicon(N_semicon,3)=Nb_pts_z_semicon;
    
    [e_semicon,z_semicon,wzs,o_semicon]=retchamp(init,struct,sh,sb,inc,{x,y},tab_semicon,[],[1:6]+7.25*i,1,1,1:6);
    [Wz,Wx,Wy]=ndgrid(wzs,wx,wy);
    W_semicon=Wz.*Wx.*Wy;
    
    E_semicon=e_semicon(:,:,:,1:3);
    H_semicon=e_semicon(:,:,:,4:6);
    x_semicon=x;y_semicon=y;
    end
    
    if trace_champ==1
        tab0=[h_air,1,Nb_pts_z+10];
        struct0={ah};
        for az=1:Nb_couches
            tab0=[tab0;[H(az),az+1,Nb_pts_z+10]];
            struct0={struct0{:},a{az}};
        end
        struct0={struct0{:},ab};
        tab0=[tab0;[h_sub,Nb_couches+2,Nb_pts_z+10]];
        tab0(tab0(:,1)>1,3)=floor(tab0(tab0(:,1)>1,1)*1000/h_2pts);
                
        [xx,wx]=retgauss(-periodicity_x/2,periodicity_x/2,15,12,[-diameter_x/2,diameter_x/2]);
        [yy,wy]=retgauss(-periodicity_y/2,periodicity_y/2,15,12,[-diameter_y/2,diameter_y/2]);
        
        if isempty(x0)==1&&isempty(z0)==1
            [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,y0},tab0,[],[1:6]+7.25*i,1,1,1:6);
        elseif isempty(y0)==1&&isempty(z0)==1
            [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{x0,yy},tab0,[],[1:6]+7.25*i,1,1,1:6);
        elseif isempty(x0)==1&&isempty(y0)==1
            tab0(:,3)=0;
            HH=cumsum(tab0(:,1));
            numz=1;
            while (HH(end)-z0)>HH(numz);numz=numz+1;end;
            tab1=[tab0(1:numz-1,:);[tab0(numz,1)-(z0-sum(tab0(numz+1:end,1))),numz,0];[0,numz,1];[z0-sum(tab0(numz+1:end,1)),numz,0];tab0(numz+1:end,:)];
            [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,yy},tab1,[],[1:6]+7.25*i,1,1,1:6);
        end
        
        for ii=1:3                                  
            o0(:,:,:,ii+3)=o0(:,:,:,ii+3)./o0(:,:,:,ii);
            o0(:,:,:,ii)=1;
        end
        indice=squeeze(sqrt(o0(:,:,:,4)));
        Ex=squeeze(e0(:,:,:,1));Ey=squeeze(e0(:,:,:,2));Ez=squeeze(e0(:,:,:,3));
        Hx=squeeze(e0(:,:,:,4));Hy=squeeze(e0(:,:,:,5));Hz=squeeze(e0(:,:,:,6));
    end
  
    %retio
    %toc(Vstart)
    %t(zou)=toc(Vstart);
    %e=mean(t);
    %pour=(zou*100./length(wavelength));
    %heure = e.*(length(wavelength)-zou)./3600;
    %minute = (heure-floor(heure))*60;
    %seconde = (minute-floor(minute))*60;
    %gilles=sprintf('%d %s  complete. Il reste environ %d heure(s), %d minute(s) et %d seconde(s)', floor(pour),'%',floor(heure), floor(minute),floor(seconde));
    %disp(gilles)
%     parfor_progress;
end
% parfor_progress(0)
 %%%% Sauvegarde des data
    %text='Si1Dpoyn_ssgranet_mm20';
    text=['diam',int2str(diameter_x*1000),'_per',int2str(periodicity_x*1000),'_h1_',int2str(h1*1000),'_h4_',int2str(h4*1000),'_h5_',int2str(h5*1000),'_Mx',int2str(Mx),'_GaAs_ISE'];
    save(text,'wavelength','ref_TE_TM','ref_TM_TM','ref_TM_TE','R0_TE_TE','R0_TE_TM','R0_TM_TM','R0_TM_TE','Abs_plots','Abs','A_sub','A_tot','test','Nb_couches','N_semicon','nh','n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11','n12','nsub','n1m','n2m','n3m','n4m','n5m','n6m','n7m','n8m','n9m','n10m','n11m','n12m','periodicity_x','periodicity_y','diameter_x','diameter_y','h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','Mx','My','pol',...
        'op_granet','Nb_pts_z','Nb_pts_z_semicon','x0','y0','z0','h_air','h_sub','Einc','Hinc','E_semicon','H_semicon','x_semicon','y_semicon','z_semicon','W_semicon','Ex','Ey','Ez','Hx','Hy','Hz','xx','yy','zz','indice','cal_abs','cal_champ','trace_champ')


%%%% Exemple de figure pour tracer le champ pour une calcul fait avec trace_champ=1, x0=[], y0=0, z0=[]
%%%% Hy en fonction de x et z
%%%% Les contours de l'objet sont en blancs (indice contient la carte des indices) 
% figure
% pcolor(xx,zz,abs(Hy).^2);shading flat;colorbar;colormap(hot);hold on;
% contour(xx,zz,indice,'w','linewidth',1.5)



