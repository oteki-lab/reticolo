function [w,disc,ddisc,z,ww]=retautomatique(d,ub,unv,z,n90,parm,nn,pmlz);
%  function [w,disc]=retautomatique(d (ou init) ,ub,unv,z,n90,parm,nn);
%
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	% construction de  maillages de tronçon   %
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  
%%%%%%%%%%%
% en 2 D  %
%%%%%%%%%%%
%  function w=retautomatique(d,ub,unv,z,n90,parm,nn);
%  d=[dx,dy]  dx dy:periode en x et y 
%  ub: texture de base (si on veut remplacer partout on peut prendre ub=[])
%
%  unv maillage de la nouvelle texture obtenu par retu dans une maille d(1) z(1) 
%
%  si n90  0 : unv est fonction de x z (option par defaut)
%  si n90  1:  unv est fonction de y z  ub est tourne de pi/2 avant
%  si n90 -1:  unv est fonction de y z  ub est tourne de- pi/2 avant
%  si n90  2:  unv est fonction de x z  ub est tourne de  pi avant
%
%  z(1) pas en z du ' reseau  fictif'  z(2)  z(3):limites en z (par defaut 0 z(1))
%  ub sera remplace par la nouvelle texture partout ou eps a des valeurs nulles
%
% si parm est indique (et eventuellement nn) trace de l'objet par rettestobjet
%
%   retautomatique construit un  MAILLAGE DE TRONCON :
%   cell array w qui permet de traiter le probleme 2D 
%   dont unv decrit une ' vue de dessus' et ub  la texture de base en xy
%    w={  { h1 ,u1}   { h2 ,u2}    ...   { hn ,un} }
%     chaque terme permettant de calculer une matrice S=retc(a,h)
%           ou a=retcouche(init,u) 
%  le premier terme est le plus haut 
%  les w peuvent etre empiles:w=[whaut,...,wmilieu,..,wbas];
%
%   il faut donc  calculer retss(S1,S2   Sn) dans cet ordre
%
%   disc{1}:discontinuites en x   (utiles a connaitre pour integrer avec gauss ...)
%   disc{2}:discontinuites en y 
%
%
%%%%%%%%%%%
% en 1 D  %
%%%%%%%%%%%
%  function [w,disc]=retautomatique(d,pol,unv,y,parm,nn);
%   disc:discontinuites en x 
%  y(1) pas en y du ' reseau  fictif'  y(2)  y(3):limites en y (par defaut 0 y(1))
%  parm nn parametres pour le trace du maillage de texture  effectivement utilisé (voir rettestobjet) 
%
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	% 'mixage' de maillages de tronçon  %
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function ww=retautomatique(w,masque);
% masque maillage d'une texture formee de n indices 1,2,  n
% w{1} maillage de tronçon a mettre dans masque = 1
%  ......
% w{n} maillage de tronçon a mettre dans masque = n
%
% si les tronçons n'ont pas la meme hauteur,ils sont prolonges vers le bas par la derniere texture
%
% ww:maillage final
%
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	% 'installation 'de pml sur un maillage de tronçon (ou de texture ) %
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  w=retautomatique(d (ou init),w,x,pmlx,y,pmly,  z,pmlz  ); en 2 D (forme statique)
%  w=retautomatique(d (ou init),w,x,pmlx,  y,pmly );        en 1 D (forme statique)
%
%  z ,pmlz en 2D et y, pmly en 1 D sont facultatifs
%  en 2D  z: cotes des limites des pml croissant.
%      l'origine est le bas,La derniere valeur de z est obligatoirement le haut du troncon
% 	les pmlz sont les pml au dessous des z (meme dimension que z)
% 	Si w est un maillage de texture, z et pmlz ont un seul element ( z doit obligatoirement etre 1)  
%   si on utilise pmlz : pas de forme dynamique 
%
%  en 1D idem avec y
%
% [w,d,disc]=retautomatique(d (ou init),w,X,pmlX,Y,pmlY,disc); en 2 D (forme dynamique)
% [w,d,disc]=retautomatique(d (ou init),w,X,pmlX,disc);        en 1 D (forme dynamique)
%
%  X discontinuitees de x (croissant)  ,pmlX: valeurs des pml en x a gauche de X 
%  Y discontinuitees de y (croissant)  ,pmlY: valeurs des pml en y a gauche de Y 
%
%  w:maillage à transformer(troncon ou texture) si w est un maillage de texture, ne pas mettre z en 2D ni y en 1D
%
% forme dynamique : les pml doivent être REELLES :
%  w est dans les coordonnees non transformees, et revient dans les coordonnees transformees
%  d est alors modifie. on peut aussi transformer un tableau disc de points de discontinuite (cell array en 2 D)
%        ( attention au centre de symetrie qui peut être modifie dans les nouvelles coordonnees ) 
%
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	% de_pml_isation d'un champ calculé par retchamp en presence de pml reelles      %
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  [e,x,y,z,{wx,wy,wz}]=retautomatique(d_phys,{ee;xx;yy;zz;{wxx,wyy,wzz}},X,pmlX,Y,pmlY ,   Z,pmlZ );   en 2 D 
%  [e,x,y,{wx,wy}]=retautomatique(d_phys,{ee;xx;yy;{wxx,wyy}},X,pmlX  ,,pmlY );            en 1 D 
%  e,x,  wx d_phys, X Y dans les coordonnes sans PML
%  ee,xx,..wxx dans les coordonnes avec PML
%    (en pratique on peut de_pml_iser avec les memes parametres qui permettent de pml_iser, d_phys est le pas avant pml)
% wx.. wxx poids pour methode de gauss (facultatif)
% transformation de points:  de l'espace reel a la pml
%  [xx,yy,zz,{wxx,wyy,wzz}]=retautomatique(d_phys,{'PML';x;y;z;{wx,wy,wz}},X,pmlX  ,,pmlY  ,   Z,pmlZ );   en 2 D 
%  [xx,yy,{wxx,wyy}]=retautomatique(d_phys,{'PML';x;y;{wx,wy}},X,pmlX  ,,pmlY );            en 1 D 
% transformation de points:  de la pml à l'espace reel
%  [x,y,z,w,{wx,wy,wz}]=retautomatique(d_phys,{'autre chaine de caracteres';xx;yy;zz;{wxx,wyy,wzz}},X,pmlX  ,,pmlY  ,   Z,pmlZ );   en 2 D 
%  [x,y,{wx,wy}]=retautomatique(d_phys,{'autre chaine de caracteres';xx;yy;{wxx,wyy}},X,pmlX  ,,pmlY );            en 1 D 
% il peut etre commode de creer des fonctions:
%                 en 1D 
% phys2num = @(x) (retautomatique(d_phys,{'PML';x;0},PML{:}));
% num2phys = @(x) (retautomatique(d_phys,{'';x;0},PML{:})); 
% Num2phys = @(ee,xx,yy) (retautomatique(d_phys,{ee;xx;yy},PML{:})); 
%                   en 2D
% phys2num = @(x,y) (retautomatique(d_phys,{'PML';x;y;0},PML{:}));
% num2phys = @(x,y) (retautomatique(d_phys,{'';x;y;0},PML{:}));
% Num2phys = @(ee,xx,yy,zz) (retautomatique(d_phys,{ee;xx;yy;zz},PML{:})); 
%
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	%  modification des eps mu dans un maillage de texture ou de tronçon %
% 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ww=retautomatique(w,ep,epnv);
% Transforme dans le maillage de texture (ou de tronçon) w, ep en epnv
%  (utile pour des milieux anisotropes dans les axes x y z)
% en 2D: ep=[k0*mu_x;k0*mu_y;k0*mu_z;k0*epsilon_x;k0*epsilon_y;k0*epsilon_z]
% pour un milieu homogene isotrope d'indice n: ep=ret2ep(n,k0);
% 
% en 1D: si pol=0  ep=[k0*epsilon;k0*mu_y;1/(k0*mu_x)]
%        si pol=2  ep=[k0*mu;k0*epsilon_y;1/(k0*epsilon_x)]
% pour un milieu homogene isotrope d'indice n: ep=retep(n,pol,k0); 
%
% pour les cylindres Popov:
% ep=[k0*mu_r;k0*mu_theta;k0*mu_z;k0*epsilon_r;k0*epsilon_theta;k0*epsilon_z]
% 
% %% Exemple
% k0=2*pi/1.5;d=[4,3];u=retu(d,{1,[0,0,2,2,3.45,1],k0});rettestobjet(d,u,0,[2,2,-d/2]);
% uu=retautomatique(u,ret2ep(3.45,k0),[1:6].');rettestobjet(d,uu,0,[2,2,-d/2]);
%
% 	%%%%%%%%%%%%%%%%%%%%%%%
% 	%  forme 'degeneree'  %
% 	%%%%%%%%%%%%%%%%%%%%%%%
%  forme 'degeneree'  [w,disc]=retautomatique(d,u,h);
%  creation directe d'un seul element w={{h,u}}
%
% See also: RETAUTO,RETTESTOBJET,RETEP,RET2EP,RETU,RETCHAMP

if nargin<3;w=retmixte(d,ub);return;end; % 'mixage' de maillages de tronçon

if nargin==3;
if numel(unv)>1;w=retchristophe(d,ub,unv); % Transforme dans le maillage de texture (ou de tronçon) w, ep en epnv
else;w={{unv,ub}};if nargout>1;disc=retdisc(d,ub);end; %  forme dégénérée
end;
return;
end;

if nargin<5;n90=[];end;
if nargin<7;nn=cell(1,2);end;
if iscell(d);d=d{end}.d;end;  % d=init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'installation 'de pml sur un maillage de tronçon  ou de_pml_isation de champ  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(unv);
    
if size(ub,1)>1;    %  d_pml_isation de champ 
if (length(ub)==3)|((length(ub)==4)&(iscell(ub{end})));    
%  [e,x,y,w]=retautomatique(d (ou init),{e;x;y;w},x,pmlx    y,pmly );            en 1 D 
if nargin<5;n90=[];parm=1;end;[w,disc,ddisc,z]=de_pml_ise1(d,ub,unv,z,n90,parm);
%  [e,x,y,z,w]=retautomatique(d (ou init),{e;x;y;z;w},x,pmlx,y,pmly  z,pmlz );   en 2 D 
else;if nargin<7;nn=cell(1,2);pmlz=1;end;[w,disc,ddisc,z,ww]=de_pml_ise2(d,ub,unv,z,n90,parm,nn,pmlz);
end;
return;end; 
% 'installation 'de pml sur un maillage de tronçon 
%  [w,d,disc]=retautomatique(d (ou init),w,x,pmlx,y,pmly,disc); en 2 D
%  [w,d,disc]=retautomatique(d (ou init),w,x,pmlx,disc); en 1 D

tr=1;if ~iscell(ub{1});if size(ub{1},1)==1;tr=0;ub={{1,ub}};else;tr=-1;end;end;% tr=1: tronçon  0 texture -1 de_pml_isation
if length(d)==1;w=ub;x=unv;pmlx=z;  % <------------------------ 1 D 
if nargin>5;y=n90;pmly=parm;ddisc=[];                       %                          <  pml en y 
ww=retautomatique(d,ub,x,pmlx);dddisc=retdisc(d,ww);ddisc=[ddisc,dddisc];
h=0;for ii=length(ww):-1:1;h=[h,h(end)+ww{ii}{1}];end;h=sort(retelimine([h,y]));hm=(h(1:end-1)+h(2:end))/2;
lh=length(hm);pml=zeros(1,lh);kw=zeros(1,lh);
h1=0;for ii=length(ww):-1:1;h2=h1+ww{ii}{1};kw((hm>h1)&(hm<h2))=ii;h1=h2;end;
h1=0;for ii=1:length(y);h2=y(ii);pml((hm>h1)&(hm<h2))=pmly(ii);h1=h2;end;
%w=cell(1,lh);for ii=1:lh;w{lh-ii+1}={(h(ii+1)-h(ii))/pml(ii),{ww{kw(ii)}{2}{1},diag([pml(ii),1/pml(ii),1/pml(ii)])*ww{kw(ii)}{2}{2}}};end;
w=cell(1,lh);for ii=1:lh;w{lh-ii+1}={h(ii+1)-h(ii),{ww{kw(ii)}{2}{1},diag([pml(ii),1/pml(ii),1/pml(ii)])*ww{kw(ii)}{2}{2}}};end;
ddisc=retelimine(ddisc);
if tr==0;w=w{1}{2};end;% modif 11 2009
return;end; %                                                                          <   pml en y

epx=zeros(3,length(pmlx));
for ii=1:length(pmlx);epx(:,ii)=pmlx(ii);end;

if nargout>1; % pml reelles 'dynamiques' nouveau pas;
ux=retu(d,x,epx);pente=real(1./ux{2}(1,:));ux=ux{1};
xx=0;XX=0;for ii=1:length(ux);xx=[xx,ux(ii)];XX=[XX,XX(end)+pente(ii)*(xx(end)-xx(end-1))];end;
disc=XX(end); % xx ancien(reel)  XX nv (apres transformation)
end;          %   pml reelles 'dynamiques'

u=retu(d,x,epx);
k=0;masque=u;uu=u{2};
w=cell(1,size(uu,2));
for ii=1:size(uu,2);k=k+1;masque{2}(:,ii)=k^2;
w{k}=ub;
for kk=1:size(w{k},2);sz=size(w{k}{kk}{2}{2});sz(1)=1;
w{k}{kk}{2}{2}=w{k}{kk}{2}{2}.*repmat(u{2}(:,ii),sz);
end;end;
w=retautomatique(w,masque);

if nargout>1;  % pml reelles 'dynamiques' modification des coordonnees
for ii=1:length(w);
w{ii}{2}{1}=retinterp(xx,XX,w{ii}{2}{1},'linear');
end;
if nargout>2;ddisc=n90;%transformation des points de discontinuite
ddisc=retinterp(xx,XX,mod(ddisc,d),'linear');
for ii=1:length(w);ddisc=[ddisc,w{ii}{2}{1}];end;% on ajoute les discontinuites de w
ddisc=sort(retelimine(ddisc));
end;% nargout>2;
end;% nargout>1  pml reelles 'dynamiques' modification des coordonnees

else;w=ub;x=unv;pmlx=z;y=n90;pmly=parm;       %  <------------------------ 2 D     

if nargin>6;if nargin==7;z=1;pmlz=nn;else;z=nn;end;nn=cell(1,2);ddisc=cell(1,2); %    <  pml en z 
ww=retautomatique(d,ub,x,pmlx,y,pmly);dddisc=retdisc(d,ww);ddisc{1}=[ddisc{1},dddisc{1}];ddisc{2}=[ddisc{2},dddisc{2}];
h=0;for ii=length(ww):-1:1;h=[h,h(end)+ww{ii}{1}];end;h=sort(retelimine([h,z]));hm=(h(1:end-1)+h(2:end))/2;
lh=length(hm);pml=zeros(1,lh);kw=zeros(1,lh);
h1=0;for ii=length(ww):-1:1;h2=h1+ww{ii}{1};kw((hm>h1)&(hm<h2))=ii;h1=h2;end;
h1=0;for ii=1:length(z);h2=z(ii);pml((hm>h1)&(hm<h2))=pmlz(ii);h1=h2;end;
w=cell(1,lh);
for ii=1:lh;
www=ww{kw(ii)}{2}{3};www(:,:,[1,2,4,5])=www(:,:,[1,2,4,5])*pml(ii);www(:,:,[3,6])=www(:,:,[3,6])/pml(ii);
%w{lh-ii+1}={(h(ii+1)-h(ii))/pml(ii),{ww{kw(ii)}{2}{1},ww{kw(ii)}{2}{2},www}};% nodif 2 2009
w{lh-ii+1}={(h(ii+1)-h(ii)),{ww{kw(ii)}{2}{1},ww{kw(ii)}{2}{2},www}};
end;
ddisc{1}=retelimine(ddisc{1});ddisc{2}=retelimine(ddisc{2});
if tr==0;w=w{1}{2};end;% modif 11 2009
return;
end; %                                            <   pml en z


epx=zeros(6,length(pmlx));epy=zeros(6,length(pmly));
for ii=1:length(pmlx);epx(:,ii)=ret2ep(1,pmlx(ii),1);end;
for ii=1:length(pmly);epy(:,ii)=ret2ep(1,1,pmly(ii));end;

if nargout>1;disc=d; % pml reelles 'dynamiques' nouveau pas;
ux=retu(d,x,epx,0,ones(6,1));pente=real(ux{3}(:,1,1));ux=ux{1};
xx=0;XX=0;for ii=1:length(ux);xx=[xx,ux(ii)];XX=[XX,XX(end)+pente(ii)*(xx(end)-xx(end-1))];end;
disc(1)=disc(1)*XX(end);XX=XX/XX(end); % xx ancien(reel)  XX nv (apres transformation)
uy=retu(d,0,ones(6,1),y,epy);pente=uy{3}(1,:,2);uy=uy{2};
yy=0;YY=0;for ii=1:length(uy);yy=[yy,uy(ii)];YY=[YY,YY(end)+pente(ii)*(yy(end)-yy(end-1))];end;
disc(2)=disc(2)*YY(end);YY=YY/YY(end);
end;               %  pml reelles 'dynamiques'


u=retu(d,x,epx,y,epy);
k=0;masque=u;uu=u{3};
w=cell(1,size(uu,1)*size(uu,2));
for ii=1:size(uu,1);for jj=1:size(uu,2);k=k+1;masque{3}(ii,jj,:)=k^2;
w{k}=ub;
for kk=1:size(w{k},2);sz=size(w{k}{kk}{2}{3});sz(3)=1;
w{k}{kk}{2}{3}=w{k}{kk}{2}{3}.*repmat(u{3}(ii,jj,:),sz);
end;end;end;

w=retautomatique(w,masque);
if nargout>1;  % pml reelles 'dynamiques' modification des coordonnees
for ii=1:length(w);
w{ii}{2}{1}=retinterp(xx,XX,w{ii}{2}{1},'linear');
w{ii}{2}{2}=retinterp(yy,YY,w{ii}{2}{2},'linear');
end;
if nargout>2;ddisc=nn;% transformation des points de discontinuite
ddisc{1}=mod([ddisc{1}/d(1),x],1);ddisc{2}=mod([ddisc{2},y]/d(2),1);
ddisc{1}=disc(1)*retelimine(retinterp(xx,XX,ddisc{1},'linear'));
ddisc{2}=disc(2)*retelimine(retinterp(yy,YY,ddisc{2},'linear'));
end;
end;           %    pml reelles 'dynamiques'

end;                   %      <------------------------ 1 D 2 D

if tr==0;w=w{1}{2};end;
return;end; % fin 'installation 'de pml  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
if length(d)==1; % 1 D 
pol=ub;disc=[];test=nargin>4;if test;if nargin>5;nn=parm;else nn=[];end;parm=n90;else;parm=[];end;n90=0;
else;            % 2 D
pol=-1;disc={[],[]};test=nargin>5;if nargin<7;nn=[];end;
if nargin<6;parm=0;end;
end;
if isempty(parm);parm=0;end;test=(parm~=0);% trace de l'objet apres 'decoupe'

if length(z)==1;z=[z,0,z];end;if length(z)==2;z=[z,z(1)];end;
if isempty(ub);ub=retu(d,0);end;% on remplace partout

[ub,d]=retrotation(n90,ub,d);%    rotation de ub de n90*pi/2

% 'mise en forme' de unv pour variation de z(2) a z(3);
n2=ceil(z(3)/z(1));
n1=floor(z(2)/z(1));

if (n2-n1)>1; % on repete unv n2-n1 fois
y=unv{2};uu=unv{3};    
if  all(all(unv{3}(:,end,:)==unv{3}(:,1,:)));y=y(1:end-1);uu=uu(:,1:end-1,:);end;

for ii=1:n2-n1-1;
unv{2}=[unv{2},y+ii];    
unv{3}=[unv{3},uu];   
end;
if unv{2}(end)~=n2-n1;unv{2}=[unv{2},n2-n1];unv{3}=[unv{3},unv{3}(:,1,:)];   end;
%if unv{2}(end)~=n2-n1+1;unv{2}=[unv{2},n2-n1+1];unv{3}=[unv{3},unv{3}(:,1,:)];   end;
unv{2}=unv{2}/(n2-n1);
%unv{2}=unv{2}/(n2-n1+1);
end;
z(2)=z(2)-n1*z(1);z(3)=z(3)-n1*z(1);z(1)=max(1,(n2-n1))*z(1);
%z(2)=z(2)-n1*z(1);z(3)=z(3)-n1*z(1);z(1)=max(1,(n2-n1+1))*z(1);

dnv=[d(1),z(1)];
epz=ret2ep(0);epun=ret2ep(1);epinf=ret2ep(inf);

if z(2)==z(3);z3=z(2)+10*eps;else;z3=z(3);end; % cas d'une epaisseur nulle 
unv=retu(dnv,0,epun,[z(2),z3],[epinf,epz],unv);% creation des bornes en z 
if test;rettestobjet(dnv,unv,parm,nn);axis equal;drawnow;end;
    
    
zz=unv{2}*z(1);xx=unv{1}*dnv(1);uu=unv{3};
uu=uu(:,:,[1,3,2,4,6,5]); % car les axes changent ...
z0=z(2);w={};
for ii=1:length(zz); if(z0>z(3));break;end;
h=min(zz(ii),z(3))-z0;z0=zz(ii);
if h>10*eps|(z(2)==z(3)); % si on veut un seul element de hauteur 0 (calcul de a)
x=xx(1);ep=reshape(uu(1,ii,:),6,1);
for iii=2:length(xx);
if ~all(reshape(uu(iii,ii,:),6,1)==ep(:,end));ep=[ep,reshape(uu(iii,ii,:),6,1)];x=[x,xx(iii)];
else;x(end)=xx(iii);end;
end;
if pol==-1; % 2 D
ww=retuu(ub,retu(d,x,ep,0,epun));
ww=retrotation(-n90,ww); %    rotation de ww de -n90*pi/2 pour retourner a l'etat initial

if nargout>1;dis=retdisc(d,ww);disc{1}=[disc{1},dis{1}];disc{2}=[disc{2},dis{2}];end; 
w=[{{h,ww}},w]; 
else;  % 1 D
if nargout>1;[ww,dis]=retu(d,x,ep);disc=[disc,dis];else;ww=retu(d,x,ep);end;    
% if pol==0;ep=ep([5,3,1],:);ep(3,:)=1./ep(3,:);w=[{{h,retu(d,x,ep)}},w];end; % E//a verifier
% if pol==2;ep=ep([2,6,4],:);ep(3,:)=1./ep(3,:);w=[{{h,retu(d,x,ep)}},w];end; % H//
if pol==0;ep=ep([6,2,1],:);ep(3,:)=1./ep(3,:);w=[{{h,retu(d,x,ep)}},w];end; % E//
if pol==2;ep=ep([3,5,4],:);ep(3,:)=1./ep(3,:);w=[{{h,retu(d,x,ep)}},w];end; % H//
end;
end;
end;

if nargout>1; % mise en forme de disc
if pol==-1;disc{1}=retelimine(sort(disc{1}));disc{2}=retelimine(sort(disc{2}));
else;disc=retelimine(sort(disc));
end;
end;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function www=retuu(w,ww);
%  function www=retuu(w,ww);
%  substitution de u dans uu la ou aucun terme de u n'est 0

x=w{1};y=w{2};u=w{3};mx=size(x,2);my=size(y,2);    
xx=ww{1};yy=ww{2};uu=ww{3};mmx=size(xx,2);mmy=size(yy,2);
[xxx,k,kk]=retelimine([x,xx]);x=xxx(kk(1:length(x)));xx=xxx(kk(length(x)+1:length(x)+length(xx)));mmmx=size(xxx,2);
[yyy,k,kk]=retelimine([y,yy]);y=yyy(kk(1:length(y)));yy=yyy(kk(length(y)+1:length(y)+length(yy)));mmmy=size(yyy,2);
% il est important de reinitialiser x xx y yy aux valeurs de xxx yyy pour les < et <= qui suivent  
uuu=zeros(mmmx,mmmy,6);    

x=[0,x];y=[0,y];xx=[0,xx];yy=[0,yy];
for ix=2:mmx+1;for iy=2:mmy+1;%on installe uuu=uu
iix=find(xxx>xx(ix-1)&xxx<=xx(ix));iiy=find(yyy>yy(iy-1)&yyy<=yy(iy));
for ii=1:6;uuu(iix,iiy,ii)=uu(ix-1,iy-1,ii);end;
end;end;

for ix=2:mx+1;for iy=2:my+1;% on substitue u a uu si aucune valeur de eps est 0
if all(u(ix-1,iy-1,:)~=0);
iix=find(xxx>x(ix-1)&xxx<=x(ix));iiy=find(yyy>y(iy-1)&yyy<=y(iy));
for ii=1:6;uuu(iix,iiy,ii)=u(ix-1,iy-1,ii);end;
end;
end;
end;

% eventuelles simplifications
nxxx=length(xxx);if nxxx>1;ixxx=[find(~all(all(uuu(2:nxxx,:,:)==uuu(1:nxxx-1,:,:),3),2)).',nxxx];else ixxx=1;end;
nyyy=length(yyy);if nyyy>1;iyyy=[find(~all(all(uuu(:,2:nyyy,:)==uuu(:,1:nyyy-1,:),3),1)),nyyy];else iyyy=1;end;
www={xxx(ixxx),yyy(iyyy),uuu(ixxx,iyyy,:)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function ww=retmixte(w,masque);
% function ww=retmixte(w,masque);
% 'mixage' de maillages de tronçon 
% masque maillage d'une texture formee de n indices 1,2,  n
% w{1} maillage de tronçon a mettre dans masque = 1
%  ......
% w{n} maillage de tronçon a mettre dans masque = n


z=cell(size(w));zmax=0;mw=length(w);
for ii=1:mw;m=length(w{ii});
z{ii}=zeros(m,1);for jj=1:m;z{ii}(jj)=w{ii}{jj}{1};end;
z{ii}=cumsum(z{ii});zmax=max(zmax,z{ii}(end));
end;
zz=0;for ii=1:length(w);z{ii}(end)=zmax;zz=[zz;z{ii}];end; % plafond commun a zmax
zz=retelimine(sort(zz));    % hauteurs de coupe
if length(zz)==1;zz=[zz,zz];end;% test
% construction de w
uu=cell(1,mw);
mzz=length(zz)-1;
ww=cell(1,mzz);
for ii=1:mzz;
ww{ii}=cell(1,2);
ww{ii}{1}=zz(ii+1)-zz(ii);
z0=(zz(ii+1)+zz(ii))/2;
for jj=1:mw; % recherche des w
k=find((z{jj}-z0)>=0);k=k(1);
uu{jj}=w{jj}{k}{2};
end;
ww{ii}{2}=retmixt(uu,masque);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uu=retmixt(u,masque);
% masque maillage d'une texture formee de n indices 1,2,  n
% u{1} a mettre dans masque=1
%  ......
% u{n} a mettre dans masque=n

if length(masque)==3;  %  2 D

xx=masque{1};yy=masque{2};
for ii=1:length(u);xx=[xx,u{ii}{1}];yy=[yy,u{ii}{2}];end;
xx=retelimine(sort(xx));yy=retelimine(sort(yy));
xxx=[0,xx];yyy=[0,yy];
m=round(sqrt(masque{3}(:,:,6)));xm=masque{1};ym=masque{2};
uu=zeros(length(xx),length(yy),6);
for ii=1:length(xx);for jj=1:length(yy);
x0=(xxx(ii+1)+xxx(ii))/2;y0=(yyy(jj+1)+yyy(jj))/2;       
im=find((xm-x0)>=0);im=im(1);
jm=find((ym-y0)>=0);jm=jm(1);
uuu=u{m(im,jm)};
iu=find((uuu{1}-x0)>=0);iu=iu(1);
ju=find((uuu{2}-y0)>=0);ju=ju(1);
uu(ii,jj,:)=uuu{3}(iu,ju,:);
end;end;

% eventuelles simplifications
nxx=length(xx);if nxx>1;ixx=[find(~all(all(uu(2:nxx,:,:)==uu(1:nxx-1,:,:),3),2)).',nxx];else ixx=1;end;
nyy=length(yy);if nyy>1;iyy=[find(~all(all(uu(:,2:nyy,:)==uu(:,1:nyy-1,:),3),1)),nyy];else iyy=1;end;

uu={xx(ixx),yy(iyy),uu(ixx,iyy,:)};

else;    %  1 D

xx=masque{1};
for ii=1:length(u);xx=[xx,u{ii}{1}];end;
xx=retelimine(sort(xx));
xxx=[0,xx];
m=round(sqrt(masque{2}(1,:)));xm=masque{1};
uu=zeros(3,length(xx));
for ii=1:length(xx);
x0=(xxx(ii+1)+xxx(ii))/2;     
im=find((xm-x0)>=0);im=im(1);
uuu=u{m(im)};
iu=find((uuu{1}-x0)>=0);iu=iu(1);
uu(:,ii)=uuu{2}(:,iu);
end;
% eventuelles simplifications
nxx=length(xx);if nxx>1;ixx=[find(~all(all(uu(:,2:nxx)==uu(:,1:nxx-1),1),1)).';nxx];else ixx=1;end;
uu={xx(ixx),uu(:,ixx)};
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   de_pml_isation du champ   %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e,x,y,w]=de_pml_ise1(d,champ,x,pmlx,y,pmly);       % 1D
ux=retu(d,x,repmat(pmlx,3,1));pente_x=real(1./ux{2}(1,:));ux=ux{1};
xx=0;XX=0;for ii=1:length(ux);xx=[xx,ux(ii)];XX=[XX,XX(end)+pente_x(ii)*(xx(end)-xx(end-1))];end;D=XX(end);
Y=champ{3};
if ~isempty(y);% pml en y
%if length(y)<length(pmly);y=[y,0];end;y(end)=max(Y);
pente_y=1./pmly;
yy=0;YY=0;for ii=1:length(y);yy=[yy,y(ii)];YY=[YY,YY(end)+pente_y(ii)*(yy(end)-yy(end-1))];end;
else;YY=[];yy=[];pente_y=[];
end;
e=champ{1};X=champ{2};if length(champ)==4;w=champ{4};else;w=[];end;
% xx yy avant pml , XX YY apres
%------------------------------------------------------
if ischar(e);  % uniquement transformee des points  
switch e;
case 'PML';% pml_isation (e est x,x est  y,y est w)
[e,X,pente_x]=transforme(xx,XX,X,1./pente_x,D,d);
[x,Y,pente_y]=transforme(yy,YY,Y,1./pente_y);
otherwise;% de_pml_isation (e est x,x est  y,y est w)
[e,X,pente_x]=transforme(XX,xx,X,pente_x,d,D);[x,Y,pente_y]=transforme(YY,yy,Y,pente_y);
end;
if ~isempty(w);w{1}=w{1}./pente_x;w{2}=w{2}./pente_y;y=w;end;
return;
end;
%------------------------------------------------------
if iscell(e);% cas des champs mis en cell array pour gagner de la memoire (imag(caly)=1 dans retchamp)
y_prv=cell(size(e));if isempty(w);w=cell(size(e));end;
for ii=1:length(e);
[e{ii},x_prv,y_prv{ii},w{ii}]=de_pml_ise1(d,{retio(e{ii});X;champ{3}{ii};w{ii}},x,pmlx,y,pmly);
e{ii}=retio(e{ii},1);
end;
x=x_prv;y=y_prv;
return;end;
%------------------------------------------------------
[x,X,pente]=transforme(XX,xx,X,pente_x,d,D);if ~isempty(w);w{1}=w{1}./pente;end;
if ~isempty(e);for ii=1:length(XX)-1;f=find((X>=XX(ii))&(X<XX(ii+1)));e(:,f,2)=e(:,f,2)*pente_x(ii);end;end;
if ~isempty(y);% pml en y
[y,Y,pente]=transforme(YY,yy,Y,pente_y);if ~isempty(w);w{2}=w{2}./pente;end;
if ~isempty(e);for ii=1:length(YY)-1;f=find((Y>=YY(ii))&(Y<YY(ii+1)));e(f,:,3:end)=e(f,:,3:end)*pente_y(ii);end;end;
else;y=champ{3};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e,x,y,z,w]=de_pml_ise2(d,champ,x,pmlx,y,pmly,z,pmlz); % 2D
D=zeros(1,2);    
ux=retu(d(1),x,repmat(pmlx,3,1));pente_x=real(1./ux{2}(1,:));ux=ux{1};
xx=0;XX=0;for ii=1:length(ux);xx=[xx,ux(ii)];XX=[XX,XX(end)+pente_x(ii)*(xx(end)-xx(end-1))];end;D(1)=XX(end);
uy=retu(d(2),y,repmat(pmly,3,1));pente_y=real(1./uy{2}(1,:));uy=uy{1};
yy=0;YY=0;for ii=1:length(uy);yy=[yy,uy(ii)];YY=[YY,YY(end)+pente_y(ii)*(yy(end)-yy(end-1))];end;D(2)=YY(end);
Z=champ{4}; 
if ~isempty(z)&~iscell(z);% pml en z
%if length(z)<length(pmlz);z=[z,0];end;z(end)=max(Z);
pente_z=1./pmlz;
zz=0;ZZ=0;for ii=1:length(z);zz=[zz,z(ii)];ZZ=[ZZ,ZZ(end)+pente_z(ii)*(zz(end)-zz(end-1))];end;
else;ZZ=[];zz=[];pente_z=[];
end;
% xx,yy avant pml , XX,YY apres
e=champ{1};X=champ{2};Y=champ{3};if length(champ)==5;w=champ{5};else;w=[];end;
%------------------------------------------------------
if ischar(e);  % uniquement transformee des points  
switch e;
case 'PML';% pml_isation   (e est x,x est  y,y est  z,z est w)
[e,X,pente_x]=transforme(xx,XX,X,1./pente_x,D(1),d(1));[x,Y,pente_y]=transforme(yy,YY,Y,1./pente_y,D(2),d(2));[y,Z,pente_z]=transforme(zz,ZZ,Z,1./pente_z);
otherwise;% de_pml_isation (e est x,x est  y,y est  z,z est w)
[e,X,pente_x]=transforme(XX,xx,X,pente_x,d(1),D(1));[x,Y,pente_y]=transforme(YY,yy,Y,pente_y,d(2),D(2));[y,Z,pente_z]=transforme(ZZ,zz,Z,pente_z);
end;
if ~isempty(w);w{1}=w{1}./pente_x;w{2}=w{2}./pente_y;w{3}=w{3}./pente_z;z=w;end;
return;
end;
%------------------------------------------------------

if iscell(e);% cas des champs mis en cell array pour gagner de la memoire
y_prv=cell(size(e));z_prv=cell(size(e));if isempty(w);w=cell(size(e));end;
for ii=1:length(e);
[e{ii},x_prv,y_prv{ii},z_prv{ii},w{ii}]=de_pml_ise1(d,{retio(e{ii});X;champ{3}{ii};w{ii}},x,pmlx,y,pmly);
e{ii}=retio(e{ii},1);
end;
x=x_prv;y=y_prv;z=z_prv;
return;end;

[x,X,pente]=transforme(XX,xx,X,pente_x,d(1),D(1));if ~isempty(w);w{1}=w{1}./pente;end;
if ~isempty(e);for ii=1:length(XX)-1;f=find((X>=XX(ii))&(X<XX(ii+1)));e(:,f,:,[1,4])=e(:,f,:,[1,4])*pente_x(ii);end;end;
[y,Y,pente]=transforme(YY,yy,Y,pente_y,d(2),D(2));if ~isempty(w);w{2}=w{2}./pente;end;
if ~isempty(e);for ii=1:length(YY)-1;f=find((Y>=YY(ii))&(Y<YY(ii+1)));e(:,:,f,[2,5])=e(:,:,f,[2,5])*pente_y(ii);end;end;

if ~isempty(z)&~iscell(z);% pml en z
[z,Z,pente]=transforme(ZZ,zz,Z,pente_z);if ~isempty(w);w{3}(f)=w{3}./pente;end;
if ~isempty(e);for ii=1:length(ZZ)-1;f=find((Z>=ZZ(ii))&(Z<ZZ(ii+1)));e(f,:,:,[3,6])=e(f,:,:,[3,6])*pente_z(ii);end;end;
else;z=champ{4};
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,X,poids]=transforme(XX,xx,X,pente,d,D);if isempty(XX);x=X;poids=ones(size(X));return;end;
if nargin<5;x=interp1(XX,xx,X,'linear');
else;pas=floor(X/D);X=X-D*pas;x=interp1(XX,xx,X,'linear')+d*pas;end;
if nargout>2;poids=ones(size(X));for ii=1:length(XX)-1;f=find((X>=XX(ii))&(X<XX(ii+1)));poids(f)=pente(ii);end;end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w=retchristophe(w,ep,epnv);


if ~iscell(w{1});w={{0,w}};tr=0;else;tr=1;end;
for ii=1:length(w);
if length(w{ii}{2})==3;             % 2D ou cylindres popov
if isempty(w{ii}{2}{3});       % cylindres popov
f=find(sum(abs(w{ii}{2}{2}-repmat(ep,1,size(w{ii}{2}{2},2))))<100*eps);
w{ii}{2}{2}(:,f)=repmat(epnv,1,length(f));
else;% 2D
prv=reshape(w{ii}{2}{3},[],6);	
f=find(sum(abs(prv-repmat(ep.',size(prv,1),1)),2)<100*eps);
prv(f,:)=repmat(retcolonne(epnv,1),length(f),1);
w{ii}{2}{3}=reshape(prv,size(w{ii}{2}{3}));
end;
else;                               % 1D
f=find(sum(abs(w{ii}{2}{2}-repmat(ep,1,size(w{ii}{2}{2},2))))<100*eps);
w{ii}{2}{2}(:,f)=repmat(epnv,1,length(f));
end;
if ~tr;w=w{1}{2};end;
end
