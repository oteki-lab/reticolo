
function varargout=retgranet(varargin);
% 
% Transformation de G.Granet appliquée sur la dérivée et non pas sur les eps mu (2011)
% semble remplacer trés avantageusement les PML reelles de Jurek
% Pour mettre en oeuvre cette transformee de coordonnee, il suffit de remplacer dans init, le parametre cao par:
% en 1D:  cao={cao_complexe;{x_disc,a_disc,b_disc,p_divise_cao}};
% en 2D:  cao={cao_complexe;{x_disc,a_disc_x,b_disc_x,p_divise_cao_x,y_disc,a_disc_y,b_disc_y,p_divise_cao_y}};
%  (Rappel si pas de transformée complexe cao_complexe=[] ou 0 )
%     x_disc centres des transformations
% Les transformations sont decrites par leur derivée dans les coordonnees physiques 
% deux types de transformations sont prevues: 'Jurek' et 'PML virtuelles'
%
% 	real(a_disc) largeur de la transformation dans les coordonnees physiques (si un seul element, il est repeté)  
%           si imag(a_disc)=0 Jurek la largeur du pic est de l'ordre de 4 a_disc
%           si imag(a_disc)=~0 PML. La largeur de la PML est a_disc avec de chaque coté une zone de passage de largeur 4/imag(a_disc)
% 	b_disc hauteur de la transformation au centre (si un seul element, il est repeté)  
% 	p_divise_cao eventuelle division de la transformee complexe  (facultatif, 1 par defaut en 2D soit on met les 2 soit on en met aucun)
%  
%  Tous les calculs tiennent compte de cette transformation: retcouche, retchamp, rets transforment eux mêmes
%  les points des coordonnees physiques aux coordonnées numériques, le centre de symetrie et le point à l'infini (cao) sont modifiés.
%  L'utilisateur travaille toujours dans les coordonnées physiques.
%  Si on désire connaitre la transformée utilisée:
% [X,Xder]=retgranet(init,x) en 1D transforme les x physiques en X numeriques et donne Xder=dX/dx
% en 2D:[X,Y,Xder,Yder]=retgranet(init,x,y), 
%       ou  xy={x,y},[XY,XYder]=retgranet(init,xy) avec X=XY{1}, Y=XY{2}, Xder=XYder{1}, Yder=XYder{2} 
% 
%  Si on désire connaitre la transformée inverse:
%  x=retgranet(init,'num2phys',X) en 1D transforme les  X numeriques en x physiques
%      [x,wx]=retgranet(init,'num2phys',X,wX)  transforme aussi les poids
% en 2D:[x,y]=retgranet(init,'num2phys',X,Y), 
%      [x,y,wx,wy]=retgranet(init,'num2phys',X,Y,wX,wY)  transforme aussi les poids
%       ou  xy={x,y},[xy,wxy]=retgranet(init,'num2phys',XY,wXY) avec x=xy{1};y=xy{2}, wx=xy{1};y=wxy{2}, X=XY{1}, Y=XY{2}, wX=wXY{1}, wY=wXY{2}
%
% %% Exemple avec tracé
% % 1D
% d=2;x_disc=[-d/4,d/4];B=100;A=.03/B;cao={[],{x_disc,A,B}};init=retinit(d,[0,0],[],[],cao);
% X=linspace(-d/2,d/2,1000);x=retgranet(init,'num2phys',X);[X,Xder]=retgranet(init,x);figure;
% subplot(2,1,1);plot(X,x,'linewidth',2);axis tight;ylabel('x physique');xlabel('X numerique');
% subplot(2,1,2);plot(X,Xder,'linewidth',2);axis tight;ylabel('dX/dx');xlabel('X numerique');retfont(gcf,0);
% 
% % 2D
% d=[1.5,3];x_disc=0;y_disc=[-d(2)/4,d(2)/4];B=100;A=.03/B;cao={[],{x_disc,A,B,y_disc,A,B}};init=retinit(d,[0,0,0,0],[],[],cao);
% X=linspace(-d(1)/2,d(1)/2,1000);Y=linspace(-d(2)/2,d(2)/2,1000);[x,y]=retgranet(init,'num2phys',X,Y);[X,Y,Xder,Yder]=retgranet(init,x,y);figure;
% subplot(2,2,1);plot(X,x,'linewidth',2);axis tight;ylabel('x physique');
% subplot(2,2,2);plot(Y,y,'linewidth',2);axis tight;ylabel('y physique');
% subplot(2,2,3);plot(X,Xder,'linewidth',2);axis tight;ylabel('dX/dx');xlabel('X numerique');
% subplot(2,2,4);plot(Y,Yder,'linewidth',2);axis tight;ylabel('dY/dy');xlabel('Y numerique');retfont(gcf,0);
%
% See also RETPASSAGE RETJUREK



% granet={x_disc,a_disc,b_disc,xmax_disc,d,d_num,X,x,wX,Xd};

if ischar(varargin{2});
switch varargin{1}{end}.dim;
case 1;[varargout{1:nargout}]=transformation_inverse_1D(varargin{:});
case 2;[varargout{1:nargout}]=transformation_inverse_2D(varargin{:});
end;	
else;
switch nargin;
case 2;[varargout{1:nargout}]=transformation_1D(varargin{:});
case 3;[varargout{1:nargout}]=transformation_2D(varargin{:});
otherwise;% initialisation
if nargout==4;[varargout{1:nargout}]=retgranet_init_1D(varargin{:});else;[varargout{1:nargout}]=retgranet_init_2D(varargin{:});end;
end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,wx,wy]=transformation_inverse_2D(init,prv,X,Y,wX,wY);
if init{end}.granet==0;x=X;y=Y;if nargin>4;wx=wX;end;if nargin>5;wy=wY;end;return;end; 
granet=init{12};
x=num2phys(X,granet{1});
if nargin>3;y=num2phys(Y,granet{2});end;
if nargin>4;[prv,prv,Xd,Yd]=retgranet(init,x,y);wx=wX./Xd;end;

if nargin>5;wy=wY./Yd;end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,wx]=transformation_inverse_1D(init,prv,X,wX);
if iscell(X);
if nargin>3;[x,y,wx,wy]=transformation_2D(init,X{:},wX{:});x={x,y};wx={wx,wy};else;[x,y]=transformation_2D(init,X{:});x={x,y};end;
return;end;
if init{end}.granet==0;x=X;if nargin>3;wx=wX;end;return;end; 
granet=init{6};
x=num2phys(X,granet);if nargin>3;[prv,Xd]=retgranet(init,x);wx=wX./Xd;end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=num2phys(X,granet);
% granet={x_disc,a_disc,b_disc,xmax_disc,d,d_num,X,x,wX,Xd};
d=granet{5};
n=floor(X/d);
x=interp1(granet{7},granet{8},X-n*d,'pchip')+n*d;
%x=interp1(granet{7},granet{8},X-n*d,'linear')+n*d;
for jj=1:10;% methode de Newton pour la transformee inverse 
prv=X-phys2num(x,granet{:});if (max(abs(prv))<1.e-10);break;end;
x=x+prv./Xderive_2x(x,granet{:});
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y,Xder,Yder]=transformation_2D(init,x,y);% transformation de coordonnees phys2num
if init{end}.granet==0;X=x;Y=y;Xder=ones(size(x));Yder=ones(size(y));return;end; 
granet=init{12};
X=phys2num(x,granet{1}{:});Y=phys2num(y,granet{2}{:});
if nargout>2;Xder=Xderive_2x(x,granet{1}{:});Yder=Xderive_2x(y,granet{2}{:});end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Xder]=transformation_1D(init,x);% transformation de coordonnees phys2num
if iscell(x);
if nargout==2;[X,Y,Xder,Yder]=transformation_2D(init,x{:});X={X,Y};Xder={Xder,Yder};else;[X,Y]=transformation_2D(init,x{:});X={X,Y};end;
return;end;
if init{end}.granet==0;X=x;Xder=ones(size(x));return;end; 
granet=init{6};
X=phys2num(x,granet{:});
if nargout>1;Xder=Xderive_2x(x,granet{:});end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cao_granet_x,cao_granet_y,granet,sym,cao]=retgranet_init_2D(x_disc,a_disc_x,b_disc_x,p_divise_cao_x,y_disc,a_disc_y,b_disc_y,p_divise_cao_y,sym,cao,d,nx,ny);
if length(sym)>1;sym_x=sym([1,3]);sym_y=sym([2,4]);else;sym_x=[];sym_y=[];end;
if length(cao)>1;cao_x=cao([1,3]);cao_y=cao([2,4]);else;cao_x=[];cao_y=[];end;
[cao_granet_x,granet_x,sym_x,cao_x]=retgranet_init_1D(x_disc,a_disc_x,b_disc_x,p_divise_cao_x,sym_x,cao_x,d(1),nx);
[cao_granet_y,granet_y,sym_y,cao_y]=retgranet_init_1D(y_disc,a_disc_y,b_disc_y,p_divise_cao_y,sym_y,cao_y,d(2),ny);
granet={granet_x,granet_y};
if ~isempty(sym_x);sym=retcolonne([sym_x;sym_y],1);end;
if ~isempty(cao_x)&~isempty(cao_y);cao=retcolonne([cao_x;cao_y],1);end;
if  isempty(cao_x)&~isempty(cao_y);cao=retcolonne([0,0;cao_y],1);end;% jamais possible
if  ~isempty(cao_x)&isempty(cao_y);cao=retcolonne([cao_x;0,0],1);end;% jamais possible
if   isempty(cao_x)&isempty(cao_y);cao=[];end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cao_granet,granet,sym,cao]=retgranet_init_1D(x_disc,a_disc,b_disc,p_divise_cao,sym,cao,d,n);
if length(a_disc)==1;a_disc=a_disc*ones(size(x_disc));end;
if length(b_disc)==1;b_disc=b_disc*ones(size(x_disc));end;
x_disc=mod(x_disc(:).',d);a_disc=a_disc(:).';b_disc=b_disc(:).';

if isempty(x_disc) | all(b_disc==0);cao_granet=zeros(1,2*n-1);cao_granet(n)=1;d_num=d;
deg=5;100;nrep=ceil(10*n/deg);[X,wX]=retgauss(0,d,deg,-nrep);granet={[],[],[],[],d,d_num,X,X,wX,ones(size(X))};return;end;% pas de transformation
f=find(b_disc~=0);x_disc=x_disc(f);a_disc=a_disc(f);b_disc=b_disc(f);

if length(cao)>1;% transformée complexe
xc=mod(real(cao(1))/p_divise_cao,d);lc=real(cao(2))/p_divise_cao;
x_disc_cao=mod([xc+lc/2+(d/p_divise_cao)*(0:p_divise_cao-1),xc-lc/2+(d/p_divise_cao)*(0:p_divise_cao-1)],d);
xmax_disc=min([abs(xc-x_disc);abs(xc+d-x_disc)],[],1)-lc/2;% zone à reserver
xmax_disc=min([abs(xc-x_disc);abs(xc+d-x_disc)],[],1)-lc/2;% zone à reserver
for ii=2:p_divise_cao;
xmax_disc=min(xmax_disc,min([abs(xc+(ii-1)*d/p_divise_cao-x_disc);abs(xc+(ii-1)*d/p_divise_cao+d-x_disc)],[],1)-1.1*lc/2);% zone à reserver
end;
else;
x_disc_cao=[];
xmax_disc=(d/2)*ones(size(x_disc));
end;
xmax_disc=(d/2)*ones(size(x_disc));% test

d_num=phys2num(d,x_disc,a_disc,b_disc,xmax_disc,d);granet={x_disc,a_disc,b_disc,xmax_disc,d,d_num};
if length(cao)>1;% modification de l'origine de la transformée complexe
cao(1)=mod(phys2num(xc*p_divise_cao,granet{:}),d);
%cao(2)=real(cao(2))*d/d_num+1i*imag(cao(2));% modif 12 2011

cao(2)=(phys2num(xc*p_divise_cao+lc/2,granet{:})-phys2num(xc*p_divise_cao-lc/2,granet{:}))*p_divise_cao+1i*imag(cao(2));

end;
if length(sym)>1;sym(2)=phys2num(sym(2),granet{:});end;% modification du centre de symetrie

% calcul de cao_granet
%deg=5;100;nrep=ceil(20*n/deg);
deg=10;nrep=ceil(10*n/deg);
[fJ,fPML]=retfind(imag(a_disc)==0);
disc=[x_disc_cao,x_disc+xmax_disc,x_disc-xmax_disc];

%for jj=[-16,-4,-1,0,1,4,16];% a tester
for jj=[-32,-8,-2,0,2,8,32];
disc=[disc,x_disc(fJ)+jj*a_disc(fJ)];
end;
for jj=[-6,-1,1,6];
disc=[disc,real(x_disc(fPML))+real(a_disc(fPML))/2+jj./imag(a_disc(fPML)),real(x_disc(fPML))-real(a_disc(fPML))/2+jj./imag(a_disc(fPML))];
end;
    
%deg=50;nrep=ceil(20*n/deg);% testif n==1;deg=1;nrep=1;end;% conique
%[X,wX]=retgauss(0,d,deg,-nrep,phys2num([x_disc_cao,x_disc,x_disc+xmax_disc,x_disc+.2*xmax_disc,x_disc-.2*xmax_disc,x_disc-xmax_disc],granet{:}),d);
%[X,wX]=retgauss(0,d,deg,-nrep,phys2num([x_disc_cao,x_disc,x_disc+xmax_disc,x_disc+.2*xmax_disc,x_disc-.2*xmax_disc,x_disc-xmax_disc],granet{:}),d);
[X,wX]=retgauss(0,d,deg,-nrep,phys2num(disc,granet{:}),d);
x_granet=linspace(0,d,2000);X_granet=phys2num(x_granet,granet{:});X_granet(end)=d;% variables provisoires  X_granet(end)=d; modit 5 2013
%for jj=1:10;% methode de Newton pour la transformee inverse 
for jj=1:10;% methode de Newton pour la transformee inverse 
x=interp1(X_granet,x_granet,X,'pchip');
prv=X-phys2num(x,granet{:});if (max(abs(prv))<1.e-10);break;end;
x=x+prv./Xderive_2x(x,granet{:});
x_granet=[x_granet,x];X_granet=[X_granet,phys2num(x,granet{:})];
[X_granet,iii]=retelimine([X_granet],-1);x_granet=x_granet(iii);
end;
x=interp1(X_granet,x_granet,X,'pchip');
%clear X_granet x_granet iii prv
Xd=Xderive_2x(x,granet{:});
if length(cao)>1;Xd_cao=ones(size(Xd));
a=imag(cao(2))/(1i-imag(cao(2)));
f=find(abs(mod(X*p_divise_cao-cao(1)+d/2,d)-d/2)<real(cao(2))/2);
prv=pi*((mod(X(f)*p_divise_cao-cao(1)+d/2,d))-d/2)/real(cao(2));
Xd_cao(f)=sin(prv).^2.*(1+a*cos(prv).^2);
else;
Xd_cao=1;
end;
cao_granet=zeros(1,2*n-1);prv=retdiag(wX.*Xd_cao.*Xd/d);
cao_granet(n:2*n-1)=sum(prv*exp((X.'*((-2i*pi/d)*[0:n-1]))),1);% pour gain de memoire
cao_granet(1:n-1)=sum(prv*exp((X.'*((-2i*pi/d)*[-n+1:-1]))),1);
granet=[granet,{X,x,wX,Xd}];

% test
% test=retcompare(exp(2i*pi/d*X.'*[-n+1:n-1])*cao_granet.',Xd_cao.*Xd)


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y=fonc1(x,r);y=(1-r)/(1+r)*(x-angle(1-r*exp(i*x)));
% function y=der_fonc1(x,r);y=(1-r)^2./(1+r^2-2*r*cos(x));
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y=fonc0(x);y=sign(x).*log(1+abs(x));
% function y=der_fonc0(x);y=1./(1+abs(x));
% function y=der2_fonc0(x);y=sign(x)./(1+abs(x)).^2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=fonc0(x,x0);if x0==0;y=asinh(x);
else;y=(((1+2*x).*atan(x0*(1+2*x))-.5*log(1+x0^2*(1+2*x).^2)/x0)-((1-2*x).*atan(x0*(1-2*x))-.5*log(1+x0^2*(1-2*x).^2)/x0))*(.5/pi);end;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=der_fonc0(x,x0);if x0==0;y=1./sqrt(1+x.^2);
else;y=(atan(x0*(1+2*x))+atan(x0*(1-2*x)))/pi;end;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=der2_fonc0(x,x0);if x0==0;y=-x./((1+x.^2).*sqrt(1+x.^2));
else;y=(1./(1+(x0*(1+2*x)).^2)-1./(1+(x0*(1-2*x)).^2))*(2*x0/pi);end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=fonc(x,xmax,x0);
b=der2_fonc0(xmax,x0)/(2*xmax);a=der_fonc0(xmax,x0)-b*xmax^2;c=der_fonc0(0,x0);
y=-a*x-(b/3)*x.^3;% partie retranchée
[f,ff]=retfind(abs(x)<xmax);
y(f)=y(f)+fonc0(x(f),x0);
[fp,fm]=retfind(x>0,ff);
prv=fonc0(xmax,x0)-a*xmax-(b/3)*xmax.^3;
y(fp)=prv;y(fm)=-prv;
y=y*(1/(c-a));% renormalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=der_fonc(x,xmax,x0);
b=der2_fonc0(xmax,x0)/(2*xmax);a=der_fonc0(xmax,x0)-b*xmax^2;c=der_fonc0(0,x0);
y=zeros(size(x));
x=abs(x);f=find(x<xmax);
y(f)=der_fonc0(x(f),x0)-a-b*x(f).^2;
y=y*(1/(c-a));% renormalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xder=Xderive_2x(x,x_disc,a_disc,b_disc,xmax_disc,d,d_num,varargin);
Xder=ones(size(x));
for ii=1:length(x_disc);
Xder=Xder+b_disc(ii)*(der_fonc((mod(x-x_disc(ii)+d/2,d)-d/2)/real(a_disc(ii)),xmax_disc(ii)/real(a_disc(ii)),imag(a_disc(ii))));
end;
Xder=Xder*d/d_num;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=inter_periodique(x,d,xmax,x0)
xx=mod(x+d/2,d)-d/2;
y=fonc(xx,xmax,x0)+2*fonc(d/2,xmax,x0)*round((x-xx)/d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X=phys2num(x,x_disc,a_disc,b_disc,xmax_disc,d,d_num,varargin);
X=x; 
for ii=1:length(x_disc);
X=X+b_disc(ii)*real(a_disc(ii))*(inter_periodique((x-x_disc(ii))/real(a_disc(ii)),d/real(a_disc(ii)),xmax_disc(ii)/real(a_disc(ii)),imag(a_disc(ii)))-inter_periodique(-x_disc(ii)/real(a_disc(ii)),d/real(a_disc(ii)),xmax_disc(ii)/real(a_disc(ii)),imag(a_disc(ii))));
end;
if nargin>6;X=X*d/d_num;end;
% % 
    