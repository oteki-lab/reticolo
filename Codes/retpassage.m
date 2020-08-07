function s=retpassage(init,varargin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construction d'une matrice s de passage granet pasgranet (sens=1) ou pasgranet granet (sens=-1) %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s=retpassage(init,ah,ab,sens);
% granet_2_pasgranet=retpassage(init_granet,ah_pasgranet,ab_granet,1);
% pasgranet_2_granet=retpassage(init_granet,ah_granet,ab_pasgranet,-1);
%( en 1D une seule des 2 matrices est utilisée, ah si sens=1, ab si sens=-1)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construction d'une matrice s de passage entre deux milieux avec des PML reelles differentes.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% en haut : inith ,texture de descripteur ah , pml construites dans les coordonnees physiques avec dphys,xh,pmlxh,( yh,pmlyh,en 2D)
% en bas  : initb ,texture de descripteur ab , pml construites dans les coordonnees physiques avec dphys,xb,pmlxb,( yb,pmlyb,en 2D)
% 1D :  s=retpassage(inith,initb,ah,ab,dphys,xh,pmlxh,xb,pmlxb);
% 2D :  s=retpassage(inith,initb,ah,ab,dphys,xh,pmlxh,yh,pmlyh,xb,pmlxb,yb,pmlyb);
% 
% Pour les reseaux avec une incidence non nulle, attention au parametre beta0 dans init
% On doit avoir beta0h*dnumh=beta0b*dnumb
%
% See also RETGRANET


if nargin>4;% PML reelles
if init{end}.dim==1;s=retpassage_1D(init,varargin{:});
else;s=retpassage_2D(init,varargin{:});
end;
else; % Granet
if init{end}.dim==1;s=retpassage_granet_1D(init,varargin{:});
else;s=retpassage_granet_2D(init,varargin{:});
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***************************%
%*      PML REELLES        *%
%***************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%  1D  %%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s=retpassage_1D(inith,initb,ah,ab,dphys,xh,pmlxh,xb,pmlxb);
ah=retio(ah);ab=retio(ab);
uh=ah{8};dh=inith{4};a1h=ah{7};betah=inith{2};sih=inith{3};
ub=ab{8};db=initb{4};a1b=ab{7};betab=initb{2};sib=initb{3};

[xxb,XXb]=caltrans_1D(dphys,xb,pmlxb);% physique  -->bas
[xxh,XXh]=caltrans_1D(dphys,xh,pmlxh);% physique  -->haut

% Points de discontinuite en haut 
XXXh=retelimine([retinterp(xxh,XXh,retinterp(XXb,xxb,ub{1},'linear'),'linear'),XXh,retinterp(xxh,XXh,xxb,'linear'),uh{1}]);
XXXb=retinterp(xxb,XXb,retinterp(XXh,xxh,XXXh,'linear'));% et leur image en bas

n=inith{1};
mub=rettestobjet(db,ab,-1,0,(XXXb(1:end-1)+XXXb(2:end))/2,3);
muh=rettestobjet(dh,ah,-1,0,(XXXh(1:end-1)+XXXh(2:end))/2,3);

Ae=calA(betab,betah,XXXb,XXXh,ones(size(mub)),1);if ~isempty(sih);Ae=sih{2}*Ae*sib{1};end;
Ah=calA(betah,betab,XXXh,XXXb,1./muh,-1)*a1h;if ~isempty(sib);Ah=sib{2}*Ah;end;

if inith{end}.sog;s={blkdiag(Ae,Ah);n;n;1};else;s={blkdiag(Ae,eye(n));blkdiag(eye(n),Ah);n;n;n;n};end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=calA(b,B,x,X,f,sens);
n=length(b);D=X(end);
A=zeros(n,n);
for ii=1:length(x)-1;
x1=x(ii);x2=x(ii+1);X1=X(ii);X2=X(ii+1);	
aa=(x2-x1)/(X2-X1);bb=(x1*X2-x2*X1)/(X2-X1);
C=repmat(-B(:),1,n)+repmat(aa*b(:).',n,1);
if sens==1;ff=f(ii);else;ff=f(ii)*aa;end;
A=A+(ff*(X2-X1))*(exp((-.5i*(X1+X2))*B(:))*exp((.5i*(X1+X2)*aa+i*bb)*b(:).')).*retsinc(C*(X2-X1)/(2*pi));
end;
A=A/D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,XX]=caltrans_1D(d,x,pmlx); 
epx=zeros(3,length(pmlx));
for ii=1:length(pmlx);epx(:,ii)=pmlx(ii);end;
ux=retu(d,x,epx);pente=real(1./ux{2}(1,:));ux=ux{1};
xx=0;XX=0;for ii=1:length(ux);xx=[xx,ux(ii)];XX=[XX,XX(end)+pente(ii)*(xx(end)-xx(end-1))];end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%  2D  %%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=retpassage_2D(inith,initb,ah,ab,dphys,xh,pmlxh,yh,pmlyh,xb,pmlxb,yb,pmlyb);
ah=retio(ah);ab=retio(ab);
sih=inith{8};sib=initb{8};betah=inith{2};betab=initb{2};n=inith{1};nx=inith{6};ny=inith{7};
%p=ah{1};pp=ah{2};q=ah{3};qq=ah{4};dh=ah{5};
if length(ah)>12;fexh=ah{8};fhxh=ah{9};feyh=ah{10};fhyh=ah{11};
uh=ah{12};dh=ah{13};else;[fexh,fhxh,feyh,fhyh,pas]=deal([]);uh=cell(1,3);end;
if length(ab)>12;fexb=ab{8};fhxb=ab{9};feyb=ab{10};fhyb=ab{11};
ub=ab{12};db=ab{13};else;[fexb,fhxb,feyb,fhyb,pas_pml]=deal([]);ub=cell(1,3);end;

[xxb,yyb,XXb,YYb]=caltrans_2D(dphys,xb,pmlxb,yb,pmlyb);% physique  -->bas 
[xxh,yyh,XXh,YYh]=caltrans_2D(dphys,xh,pmlxh,yh,pmlyh);% physique  -->haut 


% Points de discontinuite en haut 
XXXh=retelimine([retinterp(xxh,XXh,retinterp(XXb,xxb,ub{1},'linear'),'linear'),XXh,retinterp(xxh,XXh,xxb,'linear'),uh{1}]);
YYYh=retelimine([retinterp(yyh,YYh,retinterp(YYb,yyb,ub{2},'linear'),'linear'),YYh,retinterp(yyh,YYh,yyb,'linear'),uh{2}]);
XXXb=retinterp(xxb,XXb,retinterp(XXh,xxh,XXXh,'linear'));% et leur image en bas
YYYb=retinterp(yyb,YYb,retinterp(YYh,yyh,YYYh,'linear'));


eeh=rettestobjet(dh,ah,-1,0,{dh(1)*(XXXh(1:end-1)+XXXh(2:end))/2,dh(2)*(YYYh(1:end-1)+YYYh(2:end))/2},1:6);
eeb=rettestobjet(db,ab,-1,0,{db(1)*(XXXb(1:end-1)+XXXb(2:end))/2,db(2)*(YYYb(1:end-1)+YYYb(2:end))/2},1:6);
% attention à la multiplication par dh et db !!!


if length(fexh)<length(fexb);% pour gain de temps
AEx=calAA(nx,ny,dh,db,betah,betab,uh,XXXh,YYYh,XXXb,YYYb,fexh,eeh(:,:,4),-1);
AEy=calAA(nx,ny,dh,db,betah,betab,uh,XXXh,YYYh,XXXb,YYYb,feyh,eeh(:,:,5),1);
AE=blkdiag(AEx,AEy);
if ~isempty(sih);AE=sib{2}*AE*sih{1};end;AE=inv(AE);
else;
AEx=calAA(nx,ny,db,dh,betab,betah,ub,XXXb,YYYb,XXXh,YYYh,fexb,eeb(:,:,4),-1);
AEy=calAA(nx,ny,db,dh,betab,betah,ub,XXXb,YYYb,XXXh,YYYh,feyb,eeb(:,:,5),1);
AE=blkdiag(AEx,AEy);
if ~isempty(sih);AE=sih{2}*AE*sib{1};end;
end;


if length(fhxh)<length(fhxb);% pour gain de temps
AHx=calAA(nx,ny,dh,db,betah,betab,uh,XXXh,YYYh,XXXb,YYYb,fhxh,eeh(:,:,1),-1);
AHy=calAA(nx,ny,dh,db,betah,betab,uh,XXXh,YYYh,XXXb,YYYb,fhyh,eeh(:,:,2),1);
AH=blkdiag(AHx,AHy);
if ~isempty(sih);AH=sib{4}*AH*sih{3};end;
else;
AHx=calAA(nx,ny,db,dh,betab,betah,ub,XXXb,YYYb,XXXh,YYYh,fhxb,eeb(:,:,1),-1);
AHy=calAA(nx,ny,db,dh,betab,betah,ub,XXXb,YYYb,XXXh,YYYh,fhyb,eeb(:,:,2),1);
AH=blkdiag(AHx,AHy);
if ~isempty(sih);AH=sih{4}*AH*sib{3};end;AH=inv(AH);

end;	
if inith{end}.sog;s={blkdiag(AE,AH);n;n;1};else;s={blkdiag(AE,eye(n));blkdiag(eye(n),AH);n;n;n;n};end;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=calAA(nx,ny,d,D,b,B,u,x,y,X,Y,ff,f,sens);
n=nx*ny;
A=zeros(n,n);x=x*d(1);y=y*d(2);X=X*D(1);Y=Y*D(2);
xx=u{1}*d(1);yy=u{2}*d(2);mx=length(xx);my=length(yy);
Ay=cell(1,length(y)-1);[aay,bby,ky]=deal(zeros(1,length(y)-1));
Ax=cell(1,length(x)-1);[aax,bbx,kx]=deal(zeros(1,length(x)-1));

for ii=1:length(x)-1;x1=x(ii);x2=x(ii+1);X1=X(ii);X2=X(ii+1);aax(ii)=(x2-x1)/(X2-X1);bbx(ii)=(x1*X2-x2*X1)/(X2-X1);
kkk=find(xx>(x1+x2)/2);kx(ii)=kkk(1);
Cx=repmat(-B(1,1:nx).',1,nx)+repmat(aax(ii)*b(1,1:nx),nx,1);
Ax{ii}=(X2-X1)*(exp((-.5i*(X1+X2))*B(1,1:nx).')*exp((.5i*(X1+X2)*aax(ii)+i*bbx(ii))*b(1,1:nx))).*retsinc(Cx*(X2-X1)/(2*pi));
Ax{ii}=retmat(Ax{ii},ny);
end;
for jj=1:length(y)-1;y1=y(jj);y2=y(jj+1);Y1=Y(jj);Y2=Y(jj+1);aay(jj)=(y2-y1)/(Y2-Y1);bby(jj)=(y1*Y2-y2*Y1)/(Y2-Y1);
kkk=find(yy>(y1+y2)/2);ky(jj)=kkk(1);
Cy=repmat(-B(2,1:nx:n).',1,ny)+repmat(aay(jj)*b(2,1:nx:n),ny,1);
Ay{jj}=(Y2-Y1)*(exp((-.5i*(Y1+Y2))*B(2,1:nx:n).')*exp((.5i*(Y1+Y2)*aay(jj)+i*bby(jj))*b(2,1:nx:n))).*retsinc(Cy*(Y2-Y1)/(2*pi));
Ay{jj}=retmat(Ay{jj},-nx);
end;
[Bx,Kx,KKx]=retelimine(kx);
[By,Ky,KKy]=retelimine(ky);


if sens==1;%champs continus en x ( Ey Hy)
for ll=1:length(Bx);
AAA=zeros(n,n);
for ii=find(kx==Bx(ll));
AA=sparse(n,n);
for jj=1:length(y)-1;
if isempty(ff);AA=AA+Ay{jj}*aay(jj);else;AA=AA+Ay{jj}*aay(jj)/f(ii,jj);end;
end;% jj
AAA=AAA+(Ax{ii}*AA);
end;% ii
if isempty(ff);A=A+AAA;else;A=A+AAA*retio(ff{Bx(ll)});end;
end; %ll

else; %champs continus en y ( Ex Hx)
for ll=1:length(By);
AAA=zeros(n,n);
	
for jj=find(ky==By(ll));
AA=sparse(n,n);
for ii=1:length(x)-1;
if isempty(ff);AA=AA+Ax{ii}*aax(ii);else;AA=AA+Ax{ii}*aax(ii)/f(ii,jj);end;
end;% ii
AAA=AAA+(Ay{jj}*AA);
end;% jj
if isempty(ff);A=A+AAA;else;A=A+AAA*retio(ff{By(ll)});end;

end; %ll
	
end;%sens

A=A/prod(D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx,yy,XX,YY]=caltrans_2D(d,x,pmlx,y,pmly); 
epx=zeros(6,length(pmlx));epy=zeros(6,length(pmly));
for ii=1:length(pmlx);epx(:,ii)=ret2ep(1,pmlx(ii),1);end;
for ii=1:length(pmly);epy(:,ii)=ret2ep(1,1,pmly(ii));end;
ux=retu(d,x,epx,0,ones(6,1));pente=real(ux{3}(:,1,1));ux=ux{1};
xx=0;XX=0;for ii=1:length(ux);xx=[xx,ux(ii)];XX=[XX,XX(end)+pente(ii)*(xx(end)-xx(end-1))];end;
XX=XX/XX(end); % xx ancien(reel)  XX nv (apres transformation)
uy=retu(d,0,ones(6,1),y,epy);pente=uy{3}(1,:,2);uy=uy{2};
yy=0;YY=0;for ii=1:length(uy);yy=[yy,uy(ii)];YY=[YY,YY(end)+pente(ii)*(yy(end)-yy(end-1))];end;
YY=YY/YY(end);



%%%%%%%%%%%%%%%%%%%%%%%%
%**********************%
%*      GRANET        *%
%**********************%
%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%  1D  %%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s=retpassage_granet_1D(init,ah,ab,sens);
if sens==1;a=retio(ah);else a=retio(ab);end;
n=init{1};
if init{end}.granet==0;s=rets1(init);return;end;

bet=init{2};si=init{3};d=init{end}.d;[X,x,wX,Xd]=deal(init{6}{7:10});
wX=wX.'/d;
a1=a{7};% calcul de Bx dans l'espace numerique
if sens==1;% numerique en bas  au physique en haut
exp_bet_X=exp(1i*X.'*bet);exp_bet_x=exp(-1i*bet.'*x);
Ae=exp_bet_x*(retdiag(wX./Xd.')*exp_bet_X);
if all(imag(bet(:))==0);exp_bet_X=exp_bet_X';exp_bet_x=exp_bet_x';
else;
exp_bet_X=exp(-1i*bet.'*X);exp_bet_x=exp(1i*x.'*bet);
end;
if isempty(a1);
Ah=exp_bet_X*(retdiag(wX)*exp_bet_x);if ~isempty(si);Ah=Ah*si{1};end;
else;% calcul propre de Hx
X_disc=[0,a{8}{1}];% espace numerique
mu=rettestobjet(d,a,-1,0,(X_disc(1:end-1)+X_disc(2:end))/2,3);
Fac=zeros(size(wX));for ii=1:length(X_disc)-1;Fac(X>X_disc(ii) & X<X_disc(ii+1))=1/mu(ii);end;
Ah=exp_bet_X*(retdiag(wX.*Fac)*exp_bet_x)*a1;
end

else;% physique en bas  au numerique en haut
exp_bet_X=exp(-1i*bet.'*X);exp_bet_x=exp(1i*x.'*bet);
Ae=exp_bet_X*(retdiag(wX)*exp_bet_x);
if all(imag(bet(:))==0);exp_bet_X=exp_bet_X';exp_bet_x=exp_bet_x';
else;
exp_bet_X=exp(1i*X.'*bet);exp_bet_x=exp(-1i*bet.'*x);
end;
if isempty(a1);
Ah=exp_bet_x*(retdiag(wX./Xd.')*exp_bet_X);if ~isempty(si);Ah=Ah*si{1};end;
else% calcul propre de Hx
x_disc=[0,a{8}{1}];% espace physique
mu=rettestobjet(d,a,-1,0,(x_disc(1:end-1)+x_disc(2:end))/2,3);
Fac=zeros(size(wX));for ii=1:length(x_disc)-1;Fac(x>x_disc(ii) & x<x_disc(ii+1))=1/mu(ii);end;
Ah=exp_bet_x*(retdiag(wX.*Fac./Xd.')*exp_bet_X)*a1;
end
end;% sens ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(si);Ae=si{2}*Ae*si{1};end;
if ~isempty(si);Ah=si{2}*Ah;end;
if init{end}.sog;s={blkdiag(Ae,Ah);n;n;1};else;s={blkdiag(Ae,eye(n));blkdiag(eye(n),Ah);n;n;n;n};end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%  2D  %%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=retpassage_granet_2D(init,ah,ab,sens);
if init{end}.granet==0;s=rets1(init);return;end;

ah=retio(ah);ab=retio(ab);
si=init{8};bet=init{2};n=init{1};nx=init{6};ny=init{7};d=init{end}.d;
betx=bet(1,1:nx);bety=bet(2,1:nx:nx*ny);Ix=retcolonne(repmat((1:nx).',1,ny),1);Iy=retcolonne(repmat((1:ny),nx,1),1);% bet=[betx(Ix);bety(Iy)];

if length(ah)>12;fexh=ah{8};fhxh=ah{9};feyh=ah{10};fhyh=ah{11};
uh=ah{12};else;[fexh,fhxh,feyh,fhyh]=deal([]);uh=cell(1,3);end;
if length(ab)>12;fexb=ab{8};fhxb=ab{9};feyb=ab{10};fhyb=ab{11};
ub=ab{12};else;[fexb,fhxb,feyb,fhyb,pas]=deal([]);ub=cell(1,3);end;

xxxh=[0,d(1)*uh{1}];yyyh=[0,d(2)*uh{2}];% Points de discontinuite 
xxxb=[0,d(1)*ub{1}];yyyb=[0,d(2)*ub{2}];% Points de discontinuite 
eeh=rettestobjet(d,ah,-1,0,{(xxxh(1:end-1)+xxxh(2:end))/2,(yyyh(1:end-1)+yyyh(2:end))/2},1:6);% epsmu
eeb=rettestobjet(d,ab,-1,0,{(xxxb(1:end-1)+xxxb(2:end))/2,(yyyb(1:end-1)+yyyb(2:end))/2},1:6);% epsmu
[X,x,wX,Xd]=deal(init{12}{1}{7:10});[Y,y,wY,Yd]=deal(init{12}{2}{7:10});

if sens==1;% numerique en bas  au physique en haut
AEx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxb,yyyb,eeb(:,:,4),fexb,-1);% Ex est calculé en bas Integration dans l'espace physique
AEy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxb,yyyb,eeb(:,:,5),feyb,1); % Ey est calculé en bas Integration dans l'espace physique

AHx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  x,X,y,Y,wX,wY,          xxxh,yyyh,eeh(:,:,1),fhxh,-1);% Hx est calculé en haut Integration dans l'espace numerique
AHy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  x,X,y,Y,wX,wY,          xxxh,yyyh,eeh(:,:,2),fhyh,1); % Hy est calculé en haut Integration dans l'espace numerique

else;% physique en bas  au numerique en haut
AEx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  x,X,y,Y,wX,wY,          xxxb,yyyb,eeb(:,:,4),fexb,-1); % Ex est calculé en bas Integration dans l'espace numerique
AEy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  x,X,y,Y,wX,wY,          xxxb,yyyb,eeb(:,:,5),feyb,1);  % Ey est calculé en bas Integration dans l'espace numerique

AHx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxh,yyyh,eeh(:,:,1),fhxh,-1);% Hx est calculé en haut Integration dans l'espace physique
AHy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxh,yyyh,eeh(:,:,2),fhyh,1); % Hy est calculé en haut Integration dans l'espace physique

end;% sens ?
AE=blkdiag(AEx,AEy);
if ~isempty(si);AE=si{2}*AE*si{1};end;
AH=blkdiag(AHx,AHy);
if ~isempty(si);AH=si{4}*AH*si{3};end;

if init{end}.sog;% matrices S
s={blkdiag(AE,AH);n;n;1};
else;% matrices G
s={blkdiag(AE,eye(n));blkdiag(eye(n),AH);n;n;n;n};
end;

% 
%     else;
%     s={eye(2*n);blkdiag(AE,AH);n;n;n;n};
% 
%     if sens==1;% numerique en bas  au physique en haut
%     AEx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxb,yyyb,eeb(:,:,4),fexb,-1);% Ex est calculé en bas Integration dans l'espace physique
%     AEy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxb,yyyb,eeb(:,:,5),feyb,1); % Ey est calculé en bas Integration dans l'espace physique
%     AHx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxb,yyyb,eeb(:,:,1),fhxb,-1);% Hx est calculé en bas Integration dans l'espace physique
%     AHy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxb,yyyb,eeb(:,:,2),fhyb,1); % Hy est calculé en bas Integration dans l'espace physique
% 
%     else;% physique en bas  au numerique en haut
%     AEx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxh,yyyh,eeh(:,:,4),fexh,-1);% Hx est calculé en haut Integration dans l'espace physique
%     AEy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxh,yyyh,eeh(:,:,5),feyh,1); % Hy est calculé en haut Integration dans l'espace physique
%     AHx=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxh,yyyh,eeh(:,:,1),fhxh,-1);% Hx est calculé en haut Integration dans l'espace physique
%     AHy=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,  X,x,Y,y,wX./Xd,wY./Yd,  xxxh,yyyh,eeh(:,:,2),fhyh,1); % Hy est calculé en haut Integration dans l'espace physique
% 
% 
%     end;% sens ?
%     AE=blkdiag(AEx,AEy);
%     if ~isempty(si);AE=si{2}*AE*si{1};end;
%     AH=blkdiag(AHx,AHy);
%     if ~isempty(si);AH=si{4}*AH*si{3};end;
%     if sens==1;
%     s={blkdiag(AE,AH);eye(2*n);n;n;n;n};
%     else;
%     s={eye(2*n);blkdiag(AE,AH);n;n;n;n};
%     end
%     end;% matrices S ou G*********************************************
% 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,X,x,Y,y,wX,wY,xxx,yyy,ee,F,sens);
% passage des coefficients de fourier en X Y aux  coeefficients de fourier en x y integration en x y  F donne eps dans l'espace X Y , ee est epsilon en X Y  
A=zeros(nx*ny);
%for jj=1:length(F);F{jj}=1;end;ee(:)=1;%test


if sens==1;% integration en Y tranches en X (champs continus en X :Ey Hy)
eee=zeros(size(Y));
for ii=1:length(F);
for jj=1:length(yyy)-1;eee(Y>=yyy(jj) &  Y<=yyy(jj+1))=ee(ii,jj);end;% determination de ep  fonction de Y     
f=find(X>xxx(ii) &  X<=xxx(ii+1)); 
exp_bet_xX=(exp(-1i*betx.'*x(f)))  *  (retdiag(wX(f)/d(1))*exp(1i*(X(f).'*betx))) ;
exp_bet_yY=(exp(-1i*bety.'*y))  *  (retdiag(wY./(eee*d(2)))*exp(1i*(Y.'*bety)))  ; 
A=A+(exp_bet_xX(Ix,Ix).*exp_bet_yY(Iy,Iy))*F{ii};
end;% ii

else;% integration en X tranches en Y (champs continus en Y :Ex Hx )
    
eee=zeros(size(X));
for jj=1:length(F);
for ii=1:length(xxx)-1;eee(X>=xxx(ii) &  X<=xxx(ii+1))=ee(ii,jj);end;% determination de ep  fonction de X  
f=find(Y>yyy(jj) &  Y<=yyy(jj+1)); 
exp_bet_xX=(exp(-1i*betx.'*x))  *  (retdiag(wX./(eee*d(1)))*exp(1i*(X.'*betx)));
exp_bet_yY=(exp(-1i*bety.'*y(f)))  *  (retdiag(wY(f)/d(2))*exp(1i*(Y(f).'*bety)));
A=A+(exp_bet_xX(Ix,Ix).*exp_bet_yY(Iy,Iy))*F{jj};
end;% jj
end;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function A=calAA_granet(nx,ny,betx,bety,Ix,Iy,d,X,x,Y,y,wX,wY,xxx,yyy,ee,F,sens);
% % passage des coefficients de fourier en X Y aux  coeefficients de fourier en x y integration en x y  F donne eps dans l'espace X Y , ee est epsilon en X Y  
% A=zeros(nx*ny);
% %for jj=1:length(F);F{jj}=1;end;ee(:)=1;%test
% 
% 
% if sens==1;% integration en Y tranches en X (champs continus en X :Ey Hy)
% eee=zeros(size(Y));
% for ii=1:length(F);
% for jj=1:length(yyy)-1;eee(y>yyy(jj) &  y<=yyy(jj+1))=ee(ii,jj);end;% determination de ep  fonction de y     
% f=find(x>xxx(ii) &  x<=xxx(ii+1)); 
% exp_bet_xX=(exp(-1i*betx.'*x(f)))  *  (retdiag(wX(f)/d(1))*exp(1i*(X(f).'*betx))) ;
% exp_bet_yY=(exp(-1i*bety.'*y))  *  (retdiag(wY./(eee*d(2)))*exp(1i*(Y.'*bety)))  ; 
% A=A+(exp_bet_xX(Ix,Ix).*exp_bet_yY(Iy,Iy))*F{ii};
% end;% ii
% 
% else;% integration en X tranches en Y (champs continus en Y :Ex Hx )
%     
% eee=zeros(size(X));
% for jj=1:length(F);
% for ii=1:length(xxx)-1;eee(x>xxx(ii) &  x<=xxx(ii+1))=ee(ii,jj);end;% determination de ep  fonction de x  
% f=find(y>yyy(jj) &  y<=yyy(jj+1)); 
% exp_bet_xX=(exp(-1i*betx.'*x))  *  (retdiag(wX./(eee*d(1)))*exp(1i*(X.'*betx)));
% exp_bet_yY=(exp(-1i*bety.'*y(f)))  *  (retdiag(wY(f)/d(2))*exp(1i*(Y(f).'*bety)));
% A=A+(exp_bet_xX(Ix,Ix).*exp_bet_yY(Iy,Iy))*F{jj};
% end;% jj
% end;
% 

