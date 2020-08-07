function [alpha,v,f,ff,nn]=retjurek(pol,t,n,mm,teta)
% [alpha,v,f,ff,nn]=retjurek(pol,t,n,mm,teta)
% t:angles en ordre croissant sur au plus 2pi, n indices avant sur le cercle trigonométrique (meme longueur que t)
% mm nombre de termes de fourier (100 par defaut)
% teta angles où on veut calculer f
% variation de E en f(teta) r^(alpha)
% ff=d(f)/d(teta) 'bien' calculee 
% nn=indice(teta)
% pol 0 ou 2 polarisation
%  inf  pour le metal 'electrique' (si pol==0  E//=0  si pol==2  H//=0 sur le metal)
% -inf  pour le metal 'magnetique' (si pol==0  H//=0  si pol==2  E//=0 sur le metal)
%% EXEMPLE permettant de retrouver la table 4.2 p151 du Van Bladel
% for alpha=0:20:180
% nu_2=real(retjurek(2,[0,alpha]*pi/180,[1,sqrt(2)]));
% nu_5=real(retjurek(2,[0,alpha]*pi/180,[1,sqrt(5)]));
% nu_10=real(retjurek(2,[0,alpha]*pi/180,[1,sqrt(10)]));
% nu_38=real(retjurek(2,[0,alpha]*pi/180,[1,sqrt(38)]));
% nu_50=real(retjurek(2,[0,alpha]*pi/180,[1,sqrt(50)]));
% nu_100=real(retjurek(2,[0,alpha]*pi/180,[1,sqrt(100)]));
%% disp(rettexte(alpha,nu_2,nu_5,nu_10,nu_38,nu_50,nu_100));
% disp(rettexte(alpha+0,nu_2+0,nu_5+0,nu_10+0,nu_38+0,nu_50+0,nu_100+0));
% end;
%% EXEMPLE permettant de retrouver la table 4.3 p157 du Van Bladel
% for ep=[1,1.5,2,2.5,3,4,5,7.5,10,15,20,25,38,50,100];
% nu_m=retjurek(2,[90,180,360]*pi/180,[-inf,1,sqrt(ep)]);
% nu_e=retjurek(2,[90,180,360]*pi/180,[inf,1,sqrt(ep)]);
% disp(rettexte(ep,nu_m,nu_e));
% end
%
%
% CAS SPECIAL:génération des PML de Jurek
% [X,f]=retjurek(npml,pas); pas facultatif ( 2 par defaut)
%  X:de -.5,à .5 (longueur 2*npml) ,f de 1 à 1 (longueur 2*npml-1)
%  si npml=0 X=0,F=[];
% % EXEMPLE
% 		[X,f]=retjurek(3)
% % 		X = [-0.5   -0.25   -0.125    0.125    0.25    0.5]
% % 		f =       [1      0.5     0.25      0.5    1  ]
% %
% %  exemple d'utilisation:   PML={[... h0+hpml*X,...],[..., pml_avant, pml*f, ...]};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3;if nargin<2;t=2;end;[alpha,v]=jurek(pol,t);return;end;%generation des PML de Jurek
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



t=t(:).';n=n(:).';
if nargin<4;mm=100;end;
if all(isfinite(n));% tout dielectrique
u=-2*mm:2*mm;
K=diag(i*(-mm:mm));
mu=retf(t/(2*pi),n.^pol,u*2*pi);mu=inv(rettoeplitz(mu));
mmu=retf(t/(2*pi),n.^(-pol),u*2*pi);mmu=inv(rettoeplitz(mmu));

[v,alpha]=eig(-mmu*K*mu*K);
alpha=sqrt(diag(alpha));[prv,ii]=sort(real(alpha));%[prv,ii]=sort(abs(alpha));
k=2;
alpha=alpha(ii(k));v=v(:,ii(k));
if nargout>2;f=v.'*exp(i*(-mm:mm).'*teta);[prv,jj]=max(abs(f));F=f(jj);f=f/F;end;
if nargout>3;ff=(mu*K*v).'*exp(i*(-mm:mm).'*teta)/F;end;% 1/mu derivee de f
if nargout>4;% indice
t=[t(end)-2*pi,t];t0=t(1);t=t-t0;teta=mod(teta-t0,2*pi);
nn=zeros(size(teta));for ii=1:length(t)-1;nn((teta>=t(ii))&(teta<=t(ii+1)))=n(ii);end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else % un angle metallique
f=find(~isfinite(n));f=f(1);
if n(f)>0;pol_metal=1;k=1;else pol_metal=-1;k=2;end;

rota=t(f);t=t-rota;t=t(2:end);n=[n,n];n=n(f+1:f+length(t));
omega=2*pi-t(end);periode=2*(2*pi-omega);
t=[-fliplr(t),t];n=[fliplr(n(2:end)),n]; % symetrisation
alpha=(-2*mm:2*mm)*(2*pi/periode);
K=diag((2i*pi/periode)*(-mm:mm));
mu=retf(t(2:end)/periode,n.^pol,alpha*periode);mu=inv(rettoeplitz(mu));
mmu=retf(t(2:end)/periode,n.^(-pol),alpha*periode);mmu=inv(rettoeplitz(mmu));
a=-mmu*K*mu*K;

% 'reduction' de a par symetrie suivant le type de metal
if pol_metal>0 % anti symetrique le premier champ(E si pol=0 H si pol=2) est nul sur le metal
sym=sparse([1:mm,mm+2:2*mm+1],[1:mm,mm:-1:1],[ones(1,mm),-ones(1,mm)],2*mm+1,mm);ssym=.5*sym.';
else           % symetrique la derivee du premier champ(E si pol=0 H si pol=2) est nulle sur le metal
sym=sparse([1:mm+1,mm+1:2*mm+1],[1:mm+1,mm+1:-1:1],ones(1,2*mm+2),2*mm+1,mm+1);ssym=.5*sym.';ssym(mm+1,mm+1)=ssym(mm+1,mm+1)/2;
end;
a=ssym*a*sym;
[v,alpha]=eig(a);
alpha=sqrt(diag(alpha));[prv,ii]=sort(real(alpha));%[prv,ii]=sort(abs(alpha));

alpha=alpha(ii(k));v=sym*v(:,ii(k));
if nargout>2;f=v.'*exp((2i*pi/periode)*(-mm:mm).'*(teta+rota));[prv,jj]=max(abs(f));F=f(jj);f=f/F;end;
if nargout>3;ff=(mu*K*v).'*exp((2i*pi/periode)*(-mm:mm).'*(teta+rota))/F;end;% 1/mu derivee de f
if nargout>4;% indice
nnn=mod(teta+rota+periode/2,periode)-periode/2;nn=zeros(size(teta));
for ii=1:length(t)-1;nn((nnn>=t(ii))&(nnn<=t(ii+1)))=n(ii);end;
end;
end;% un angle metallique ?	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,f]=jurek(npml,pas);%generation des PML de Jurek
if npml==0,X=0;f=[];return;end;
X=.5./(pas.^(npml-1:-1:0));X=[-fliplr(X),X];
f=1./(pas.^(npml-1:-1:0));f=[fliplr(f),f(2:end)];
