function filename=reticolo_eng(in)
%%%%                Calculate diffraction in a system composed of n layers
%%%%                                         (22/09/2020)
%%%%
%%%% The code is originally implemented to have rectangular structures in a
%%%% rectangular lattice, where the size and the period of the structures
%%%% is the same for all layers which are structured.
%%%% The period is supposed inferior to the wavelength, such that only the
%%%% specular reflection is calculated.
%%%% The substrate is supposed to absorb all remaining light, such that we
%%%% do not calculate transmission
%%%%
%%%% Calculated for all wavelengths
%%%% - Reflection coefficient (giving the total absorption A=1-R)
%%%% - Absorption in each layer
%%%%
%%%% Can also calculate, for a given wavelength
%%%% - The volumic EM field in a given layer
%%%% - A cross section of the EM field along a given plane
%%%%
%%%% Unit of length; everything is in ??m
%%%%

% -----------------------------------------
% The layers are numbered with i : refractive indices ni et nim ; thickness hi.
% nim is the index inside the structure (see side view), ni is the index
% between the structures.
% IMPORTANT : If a given layer does not have structures, write nim=0
% The number of layers is NB_couches, they are numbered from 1 to 12,
% starting from the top one
% If there is one layer, it must be layer one, if two layers, 1 and 2, etc
% IMPORTANT : for unused layers, write hi=0
%
% The incident medium has an index nh and must be transparent (nh real)
% Le substrat est un metal (or an absorbing medium) with index nsub

% -----------------------------------------
% Stack (view from the side)
%
%          Incident medium (nh)
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
%           Metallic substrate (nsub)
% z ^
%   |
%   --->  x

% -----------------------------
% View from above
%
%                              period_x
%                             <------------->
%         *******       *******       *******       *******       ******* ^
%         *******       *******       *******       *******       ******* |
%         *******       *******       *******       *******       ******* | period_y
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
% Outputs (vectors are the same size as wavelength, matrices are of size Nb_couches x length(wavelength), except the fields)
% ------------------------------------------------------------------------------------------------------------------------------------
% Only the specular reflection is calculated (the period is supposed inferior to the wavelength, so that the higher orders are in total internal reflection)
%
% ref_TE_TE : reflection coefficient for the amplitude of the field for an incident wave TE and a reflected wave TE
% ref_TE_TM : reflection coefficient for the amplitude of the field for an incident wave TE and a reflected wave TM
% ref_TM_TM : reflection coefficient for the amplitude of the field for an incident wave TM and a reflected wave TM
% ref_TM_TE : reflection coefficient for the amplitude of the field for an incident wave TM and a reflected wave TE
% R0_TE_TE : reflection coefficient for the intensity for TE --> TE
% R0_TE_TM : reflection coefficient for the intensity for TE --> TM
% R0_TM_TM : reflection coefficient for the intensity for TM --> TM
% R0_TM_TE : reflection coefficient for the intensity for TM --> TE
%
% Note : The four reflection coefficient are different from zero only if
% both angles are different from zero
% As long as one of the angles is zero, the programme uses the symmetry of
% the problem to get only one non-zero coefficient
% The polarization is then fixed with pol
%
% If cal_abs=1
% --------------
% Abs_plots(i,:) : Absorption in the structures of the layer i, as a
% function of the wavelength
% Abs(i,:) : Absorption in the totality of the layer i, as a
% function of the wavelength
% A_sub : absorption in the substrate
% A_tot : total absorption
% test : 1-R0-A_tot should be close to zero (indicator of poor convergence)
%
% If cal_champ=1  (only one wavelength)
% ----------------------------------------------------------
% E_semicon, H_semicon : E and H fields in the layer specified by N_semicon
%                        The field is calculated in points (x_semicon,y_semicon,z_semicon)
%                        dimension of the vectors (Nb_pts_z_semicon,165,150,3)
%                        E_semicon(:,:,:,1) is Ex for all the points
% W_semicon : to calculate the integral of E_semicon (Gauss method)
%             int_{maille}|E|^2 = sum(sum(sum(W_semicon.*sum(abs(E_semicon).^2,4))))
% Einc, Hinc : E and H incident fields used to calculate E_semicon and H_semicon
%
% If trace_champ=1  (only one wavelength)
% ------------------------------------------------------------
% Ex,Ey,Ez : cross-section of the E field at points zz,xx,yy
% Hx,Hy,Hz : cross-section of the H field at points zz,xx,yy
% The direction of the cross-section can be chosen among x=x0, y=y0 ou z=z0

%[prv,vmax]=retio([],inf*1i);

%% Parameters of the structure and the calculation
% Wavelengths and angle of incidence
npoints=in.npoints;                                 % point number of wavelength
lambdamin=in.lambdamin;                             % min wavelength
lambdamax=in.lambdamax;                             % max wavelength
wavelength=linspace(lambdamin,lambdamax,npoints);   % range of wavelength
theta=[0,0];                                        % angle of incidence in degrees

% Geometric parameters
period_x=in.period_x;                               % period in x
period_y=in.period_y;                               % period in y

% Number of Fourier terms
Mx=double(in.Mx);                                   % Number of Fourier terms in x
My=double(in.My);                                   % Number of Fourier terms in y

% Parameters of each layer
nh=1;                                               % Refraction indices of Air (front)
nsub=ones(size(wavelength));                        % Refraction indices of Air (back)
n_layer = size(in.params,1);                     % Total number of layers

%% Numerical parameters
pol=in.pol;                 % For normal incidence, TM <=> H//y and TE <=> E//y

% IMPORTANT: To be changed if non-normal incident or if non-rectangular structures
sym = tif(in.sym, [pol-1,pol-1,0,0], []); % [pol-1,pol-1,0,0]: The symmetry of the structure, more symmetry means shorter calculation time
% if theta(1)==0 && theta(2)==0;sym=[pol-1,pol-1,0,0];end;
% if theta(1)~=0 && theta(2)==0;sym=[0,pol-1,0,0];end;
% if theta(1)==0 && theta(2)~=0;sym=[1-pol,0,0,0];end;
% if theta(1)~=0 && theta(2)~=0;sym=[];end;

cal_abs = npoints>1;        % If 1, calculate absorption in each layer
Nb_pts_z=10;                % Number of points in z to calculate absorption, when absorption is calculated in each layer

% IMPORTANT: only one wavelength (calculate can be quite heavy, depending on Nb_pts_z_semicon)
cal_champ=0;                % If 1, calculate the field in layer N_semicon
N_semicon=3;                % Layer where field is calculated (if cal_champ=1)
Nb_pts_z_semicon=50;        % Number of points in the z direction to calculate the field in N_semicon (if cal_champ=1)

% IMPORTANT: only one wavelength
trace_champ = in.trace_champ;%npoints==1;   % si 1,calculates a cross-section of the field
x0 = in.cs_x;               % Cross section along x=x0 if trace_champ=1 ([] if the cross-section is along another direction)
y0 = in.cs_y;               % Cross section along y=y0 if trace_champ=1 ([] if the cross-section is along another direction)
z0 = in.cs_z;               % Cross section along z=z0 if trace_champ=1 ([] if the cross-section is along another direction)
                            % note: z=0 corresponds to the bottom of the considered stack, at a depth h_sub inside the substrate
h_air=0.05;                 % Thickness in incident medium to represent the cross-section (trace_champ=1)
h_sub=0.05;                 % Thickness in the substrate to represent the cross-section  (trace_champ=1)
h_2pts=20;                  % distance between 2 points in z pour the layers whose thickness is higher than 1??m (trace_champ=1)
op_objet=0;                 % If 1, plot the geometry to verify the calculated structure is correct

% IMPORTANT: this parameter is tricky to use, and does not work out of normal incidence. Better keep it at zero
op_granet=0;                % If 1, RCWA is modified to improve convergence (Transforms the real coordinates at discontinuities)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code (for parallel computing of the absorption)

% Initialisation of the parameters for parallel computing
ref_TE_TE=zeros(1,length(wavelength)); R0_TE_TE=zeros(1,length(wavelength));
ref_TE_TM=zeros(1,length(wavelength)); R0_TE_TM=zeros(1,length(wavelength));
ref_TM_TM=zeros(1,length(wavelength)); R0_TM_TM=zeros(1,length(wavelength));
ref_TM_TE=zeros(1,length(wavelength)); R0_TM_TE=zeros(1,length(wavelength));
A_tot=zeros(1,length(wavelength));
A_sub=zeros(1,length(wavelength));
Abs=zeros(n_layer,length(wavelength));
Abs_plots=zeros(n_layer,length(wavelength));
XX=[]; YY=[]; ZZ=[]; E=[]; I=[]; CONTOUR=[]; CS=[];
Ntre=1;
H=cell2mat(in.params(:,1));
n = cell2mat(in.params(:,2));
nm = cell2mat(in.params(:,3));
active_layer = in.layers{strcmp('Active region', in.layers),2};
active_t = sum(H(min(active_layer):end));
active_b = sum(H(max(active_layer)+1:end));
if cal_abs||cal_champ==1||trace_champ||in.cal_field;op_retcouche=1; else; op_retcouche=0; end
if H(n_layer)<1e-5;disp('WARNING : There is a problem in the definition of the layers number !!'); return; end
if trace_champ&&isempty(x0)==1&&isempty(y0)==1&&isempty(z0)==1;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end
if trace_champ&&isempty(x0)==0&&isempty(y0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end
if trace_champ&&isempty(x0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end
if trace_champ&&isempty(y0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end

parfor zou=1:length(wavelength)
    disp(['Calculation n-' int2str(zou) ' of ' int2str(length(wavelength))])

    inc=[];
    e0=[];
    o0=[];
    R0_TE_TE_vect=zeros(1); ref_TE_TE_vect=zeros(1);
    R0_TE_TM_vect=zeros(1); ref_TE_TM_vect=zeros(1);
    R0_TM_TM_vect=zeros(1); ref_TM_TM_vect=zeros(1);
    R0_TM_TE_vect=zeros(1); ref_TM_TE_vect=zeros(1);
    A_tot_vect=zeros(1);
    A_sub_vect=zeros(1);
    Abs_vect=zeros(n_layer,1);
    Abs_plots_vect=zeros(n_layer,1);
    test_vect=zeros(1);
    xx=[]; yy=[]; zz=[]; Ex=[]; Ey=[]; Ez=[]; Hx=[]; Hy=[]; Hz=[]; indice=[]; cs=0;
    Einc=[]; E_semicon=[]; x_semicon=[]; z_semicon=[];
    Hinc=[]; H_semicon=[]; y_semicon=[]; W_semicon=[];
    
    lwa=wavelength(zou);
    k0=2*pi/lwa;
    period=[period_x, period_y];

    ns=nsub(zou);
    N=set_layer_number(n,zou,n_layer);      Nm=set_layer_number(nm,zou,n_layer);
    diameter_x=cell2mat(in.params(:,4));    diameter_y=cell2mat(in.params(:,5));
    xdisc=[-diameter_x/2,diameter_x/2];     ydisc=[-diameter_y/2,diameter_y/2];
    Bx=500; Ax=0.02/Bx; By=Bx; Ay=Ax;
    beta=[k0*nh*sin(theta(1)*pi/180), k0*nh*sin(theta(2)*pi/180)];

    uh=retu(period,{nh,k0});
    ub=retu(period,{ns,k0});

    if op_granet==1
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym,{[],{xdisc,Ax,Bx,ydisc,Ay,By}});
        ah=retcouche(init,uh,1);
    else
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym);
        ah=retcouche(init,uh,op_retcouche);
    end
    ab=retcouche(init,ub,op_retcouche);

    % get structure
    [u,a] = structure(in,n_layer,period,N,Nm,k0,diameter_x,diameter_y,Ntre,init,op_retcouche);

    if op_objet==1
        struct_test=cell(1,n_layer+2);
        struct_test{1}={0.1,uh};
        for az=1:n_layer
            struct_test{az+1}={H(az),u{az}};
        end
        struct_test{n_layer+2}={0.05,ub};
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
    K=beta+ordre*2*pi./period;
    inch=find(((K(1)-beth(:,1)).^2+(K(2)-beth(:,2)).^2)<1e-8);
    difh=inch;
    sh=rettronc(sh,haut(inch,1),haut(difh,1),1);

    [sb,bas,betb,cb,anglb]=retb(init,ab,-0.001);

    incb=[]; difb=[];
    sb=rettronc(sb,bas(incb,1),bas(difb,1),-1);
    stemp=retc(a{1},H(1));
    for az=2:n_layer
        stemp=retss(stemp,retc(a{az},H(az)));
    end

    stot=retss(sh,stemp,sb);
    ef=retreseau(init,stot,betb,cb,anglb,incb,difb,beth,ch,anglh,inch,difh);
    if isempty(ef.amplitude(ef.dif.TEh,ef.inc.TEh))==0;ref_TE_TE_vect=ef.amplitude(ef.dif.TEh,ef.inc.TEh); end
    if isempty(ef.amplitude(ef.dif.TMh,ef.inc.TEh))==0;ref_TE_TM_vect=ef.amplitude(ef.dif.TMh,ef.inc.TEh); end
    if isempty(ef.amplitude(ef.dif.TMh,ef.inc.TMh))==0;ref_TM_TM_vect=ef.amplitude(ef.dif.TMh,ef.inc.TMh); end
    if isempty(ef.amplitude(ef.dif.TEh,ef.inc.TMh))==0;ref_TM_TE_vect=ef.amplitude(ef.dif.TEh,ef.inc.TMh); end

    R0_TE_TE_vect=abs(ref_TE_TE_vect)^2;
    R0_TE_TM_vect=abs(ref_TE_TM_vect)^2;
    R0_TM_TM_vect=abs(ref_TM_TM_vect)^2;
    R0_TM_TE_vect=abs(ref_TM_TE_vect)^2;

    R0_TE_TE(:,zou)=R0_TE_TE_vect;
    R0_TE_TM(:,zou)=R0_TE_TM_vect;
    R0_TM_TE(:,zou)=R0_TM_TE_vect;
    R0_TM_TM(:,zou)=R0_TM_TM_vect;
    ref_TE_TE(:,zou)=ref_TE_TE_vect;
    ref_TE_TM(:,zou)=ref_TE_TM_vect;
    ref_TM_TE(:,zou)=ref_TM_TE_vect;
    ref_TM_TM(:,zou)=ref_TM_TM_vect;

    if op_granet==1
        % [Xdisc,Ydisc]=retgranet(init,[-diameter_x/2,-diameter_x/2+w_rectangle,diameter_x/2,diameter_x/2+h_rectangle],[-diameter_y/2,-diameter_y/2+w_rectangle,diameter_y/2,diameter_y/2+h_rectangle]);
        [Xdisc,Ydisc]=retgranet(init,[-diameter_x/2,diameter_x/2],[-diameter_y/2,diameter_y/2]);
        [X,wX]=retgauss(-period_x/2,period_x/2,15,10,Xdisc);
        [Y,wY]=retgauss(-period_y/2,period_y/2,15,10,Ydisc);
        [x,y]=retgranet(init,'num2phys',X,Y);
        [X,Y,Xder,Yder]=retgranet(init,x,y);
        wx=wX./Xder; wy=wY./Yder;
    else
        [x,wx]=retgauss(-period_x/2,period_x/2,15,10,[-unique(diameter_x)/2,unique(diameter_x)/2]);
        [y,wy]=retgauss(-period_y/2,period_y/2,15,10,[-unique(diameter_y)/2,unique(diameter_y)/2]);
    end

    if cal_abs||cal_champ==1||trace_champ||in.cal_field
        sb_norm=retb(init,ah,-0.1,0,[],[]);
        tab_norm=[0,1,1];
        if size(ef.inc.teta,2)==1
            inc=1;
        elseif size(ef.inc.teta,2)==2
            inc=[0,0];
            if pol==0;inc(ef.inc.TE)=1; elseif pol==2; inc(ef.inc.TM)=1; end
        end
        for az=1:(n_layer)
            [ei,zi]=retchamp(init,{ah},sh,sb_norm,inc,{x,y},tab_norm,[],1:6,1,1,1:6);
            flux=retpoynting(ei,[0,0,-1],wx,wy,[]);
            inc=1/sqrt(flux).*inc;

            [einc,zinc]=retchamp(init,{ah},sh,sb_norm,inc,{0,0},tab_norm,[],1:6,1,1,1:6);
            Einc=squeeze(einc(:,:,:,1:3));
            Hinc=squeeze(einc(:,:,:,4:6));
        end
    end

    if cal_abs
        nunu=find(Nm~=0);
        num=[];
        for zer=1:length(nunu)
            if nunu(zer)==1
                num=[num,nunu(zer),nunu(zer)+1];
            else
                num=[num,nunu(zer)-1,nunu(zer),nunu(zer)+1];
            end
        end
        num=retelimine(num);

        tab=zeros(n_layer,3);
        tab2=[];
        struct={ab};
        for az=1:n_layer
            tab(az,:)=[H(az),az+1,Nb_pts_z];               %tab(az,:)=[H(az),az+1,0];
            tab2=[tab2;[0,az+1,1];[H(az),az+1,0]];
            struct=[struct(:)',a(az)];
        end
        tab(num,3)=Nb_pts_z;
        tab2=[tab2;[0,1,1];[1,1,0];[0,1,1]];

        [e,z,wz,o]=retchamp(init,struct,sh,sb,inc,{x,y},tab,[],(1:6)+7.25i,1,1,1:6);
        [e2,z2,wz2,o2]=retchamp(init,struct,sh,sb,inc,{x,y},tab2,[],(1:6)+7.25i,1,1,1:6);

        for ii=1:3
            o(:,:,:,ii+3)=o(:,:,:,ii+3)./o(:,:,:,ii);
            o(:,:,:,ii)=1;
        end
        [Wz,Wx,Wy]=ndgrid(wz,wx,wy);
        W=Wz.*Wx.*Wy;
        flux_poyn=retpoynting(e2,[0,0,-1],wx,wy,[]);

        A_sub_vect=flux_poyn(2)-flux_poyn(1);

        for az=1:n_layer
            Abs_vect(n_layer-az+1)=flux_poyn(az+2)-flux_poyn(az+1);
        end

        for az=1:length(num)
            numx=find(x>-diameter_x(az)/2 & x<diameter_x(az)/2);
            numy=find(y>-diameter_y(az)/2 & y<diameter_y(az)/2);
            Abs_plots_vect(num(az))=0.5*k0*sum(sum(sum(W((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,numx(az),numy(az)).*imag(o((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,numx(az),numy(az),4)).*sum(abs(e((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,numx(az),numy(az),1:3)).^2,4))));
            Abs_vect(num(az))=0.5*k0*sum(sum(sum(W((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,:,:).*imag(o((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,:,:,4)).*sum(abs(e((length(num)-az)*Nb_pts_z+1:(length(num)-az+1)*Nb_pts_z,:,:,1:3)).^2,4))));
        end

        Abs(:,zou)=Abs_vect;
        Abs_plots(:,zou)=Abs_plots_vect;
        A_tot_vect=flux_poyn(end)-flux_poyn(1);
        test_vect=1-R0_TE_TE_vect-R0_TE_TM_vect-R0_TM_TM_vect-R0_TM_TE_vect-A_tot_vect;
        A_tot(:,zou)=A_tot_vect;
        A_sub(:,zou)=A_sub_vect;
    end

    if cal_champ==1
        tab_semicon=zeros(n_layer);  %tab_semicon=[];
        struct={ab};
        for az=1:n_layer
            tab_semicon(az)=[tab_semicon;[H(az),az+1,0]];   %tab_semicon=[tab_semicon;[H(az),az+1,0]];
            struct=[struct(:)',a(az)];
        end
        tab_semicon(N_semicon,3)=Nb_pts_z_semicon;

        [e_semicon,z_semicon,wzs,o_semicon]=retchamp(init,struct,sh,sb,inc,{x,y},tab_semicon,[],(1:6)+7.25i,1,1,1:6);
        [Wz,Wx,Wy]=ndgrid(wzs,wx,wy);
        W_semicon=Wz.*Wx.*Wy;

        E_semicon=e_semicon(:,:,:,1:3);
        H_semicon=e_semicon(:,:,:,4:6);
        x_semicon=x; y_semicon=y;
    end

    if trace_champ||in.cal_field
        tab0 = zeros(n_layer+2,3);
        tab0(1,:) = [h_air,1,Nb_pts_z+10];
        struct0={ah};
        for az=1:n_layer
            tab0(az+1,:)=[H(az),az+1,Nb_pts_z+10];
            struct0=[struct0(:)',a(az)];
        end
        struct0=[struct0(:)',{ab}];
        tab0(n_layer+2,:)=[h_sub,n_layer+2,Nb_pts_z+10];
        tab0(tab0(:,1)>1,3)=floor(tab0(tab0(:,1)>1,1)*1000/h_2pts);

        [xx,wx]=retgauss(-period_x/2,period_x/2,15,12,[-diameter_x/2,diameter_x/2]);
        [yy,wy]=retgauss(-period_y/2,period_y/2,15,12,[-diameter_y/2,diameter_y/2]);

        if trace_champ
            if isempty(x0)==1&&isempty(z0)==1
                cs=y0;
                [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,y0},tab0,[],(1:6)+7.25i,1,1,1:6);
            elseif isempty(y0)==1&&isempty(z0)==1
                cs=x0;
                [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{x0,yy},tab0,[],(1:6)+7.25i,1,1,1:6);
            elseif isempty(x0)==1&&isempty(y0)==1
                cs=z0;
                tab0(:,3)=0;
                HH=cumsum(tab0(:,1));
                numz=1;
                while (HH(end)-z0)>HH(numz);numz=numz+1; end
                tab1=[tab0(1:numz-1,:);[tab0(numz,1)-(z0-sum(tab0(numz+1:end,1))),numz,0];[0,numz,1];[z0-sum(tab0(numz+1:end,1)),numz,0];tab0(numz+1:end,:)];
                [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,yy},tab1,[],(1:6)+7.25i,1,1,1:6);
            end

            for ii=1:3
                o0(:,:,:,ii+3)=o0(:,:,:,ii+3)./o0(:,:,:,ii);
                o0(:,:,:,ii)=1;
            end
            indice=squeeze(sqrt(o0(:,:,:,4)));
            Ex=squeeze(e0(:,:,:,1)); Ey=squeeze(e0(:,:,:,2)); Ez=squeeze(e0(:,:,:,3));
            Hx=squeeze(e0(:,:,:,4)); Hy=squeeze(e0(:,:,:,5)); Hz=squeeze(e0(:,:,:,6));
            XX=[XX,xx]; YY=[YY,yy]; ZZ=[ZZ,zz];
            E = [E, tif(pol==0, Ey, Ex)];

            ct = repmat(cs,size(indice));
            ct(indice==1) = -1*(period_x+period_y);
            CONTOUR = [CONTOUR,ct];
            CS = [CS,cs];
        else
            tab0(:,3)=0;
            HH=cumsum(tab0(:,1));
            numz=1;

            zz = sum(H)-active_b:-0.001:sum(H)-active_t;
            for zzi = zz
                while (HH(end)-(sum(H)+h_sub-zzi))>HH(numz);numz=numz+1; end
                tab1=[tab0(1:numz-1,:);[tab0(numz,1)-((sum(H)+h_sub-zzi)-sum(tab0(numz+1:end,1))),numz,0];[0,numz,1];[(sum(H)+h_sub-zzi)-sum(tab0(numz+1:end,1)),numz,0];tab0(numz+1:end,:)];
                [e0,~,~,~]=retchamp(init,struct0,sh,sb,inc,{xx,yy},tab1,[],(1:6)+7.25i,1,1,1:6);
                Ez=cat(3,Ez,squeeze(e0(:,:,:,tif(pol==2, 1, 2))));
            end
            ZZ=[ZZ,zz];
            mapping_field(in, zou, xx, yy, zz, Ez);
            I(:,zou) = mean(abs(Ez).^2, [1 2]);
        end
    end
end

%% Saving and plotting output data
if in.cal_field
    text=append(in.prefix, 'I_mean.mat');
    save(text, 'I');
    figure
    hP=surf(wavelength, ZZ, I);
    shading interp;
    colormap(jet);
    colorbar;
    caxis([0 Inf]);
    view(2)
    hP.DataTipTemplate.DataTipRows(1).Label = 'Wavelength';
    hP.DataTipTemplate.DataTipRows(2).Label = 'Depth';
    hP.DataTipTemplate.DataTipRows(3).Label = '|E|^2';
    hP.DataTipTemplate.DataTipRows(3).Value = hP.CData;
    filename = append(in.prefix,"I_mean.png");
    saveas(gcf, filename);
end

%%%% Plot a cross section with trace_champ=1, x0=[], y0=0, z0=[]
if trace_champ
    if isempty(x0)==1&&isempty(z0)==1
        filename = save_cross_section(in, XX, ZZ, CONTOUR, CS(1), E, 'x', 'z', 'y', [-period_x/2,period_x/2], [0 inf], [-period_y/2,period_y/2]);
    elseif isempty(y0)==1&&isempty(z0)==1
        filename = save_cross_section(in, YY, ZZ, CONTOUR, CS(1), E, 'y', 'z', 'x', [-period_y/2,period_y/2], [0 inf], [-period_x/2,period_x/2]);
    elseif isempty(x0)==1&&isempty(y0)==1
        filename = save_cross_section(in, YY, XX, CONTOUR, CS(1), E, 'y', 'x', 'z', [-period_y/2,period_y/2], [-period_x/2,period_x/2], [0 inf]);
    end
end

%%%% Plot Absorption - wavelength
if cal_abs
    %%%% Save the data into a file
    text=append(in.prefix, in.res, '.mat');
    
    layers_name = in.layers(:,1);
    layer_numss = in.layers(:,2);
    legends = cell(1,length(in.layers));
    Abs_array = cell(1,length(in.layers));
    for index=1:length(in.layers)
        layer_nums = layer_numss{index};
        if length(layer_nums)>1
            legends(index) = append(layers_name(index),'(',int2str(layer_nums(1)),'-',int2str(layer_nums(length(layer_nums))),')');
        else
            legends(index) = append(layers_name(index),'(',int2str(layer_nums(1)),')');
        end
        Abs_temp = zeros(1,length(wavelength));
        for index2=1:length(layer_nums)
            Abs_temp = Abs_temp + Abs(layer_nums(index2),:);
        end
        Abs_array{index} = Abs_temp;
    end

    %%%% Example to plot the absorption
    figure
    hold on
    for index=1:length(in.layers)
        plot(wavelength,cell2mat(Abs_array(:,index)), 'Linewidth',3);
    end
    hold off
    legend(legends)
    xlabel('\lambda (um)')
    ylabel('Absorption')
    xlim([min(wavelength) max(wavelength)])
    ylim([0 1])
    set(gca,'Fontsize',12)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    box on

    filename = append(in.prefix,"absorption graph.png");
    saveas(gcf, filename);
    
    % constants
    E_G=1.424;          % bandgap
    h=6.63e-34;         % planck constant
    c=3e8;              % speed of light
    q=1.60e-19;         % elementary charge

    wavelength_A = wavelength;
    load('Solarspectrum.mat','wavelength')                                          % spectra loaded in [W.m^-2.nm^-1]
    w_range = find(wavelength==lambdamin*1000):find(wavelength==lambdamax*1000+1);
    wavelength_S = wavelength(w_range);                                             % wavelength from wavmin to wavmax [nm]
    E_A = fliplr(h*c./(q*wavelength_A*1e-6));                                       % from wavelength [nm] to energy [eV]
    E_S = fliplr(flipud(h*c./(wavelength_S*1e-9))/q);                               % from wavelength [nm] to energy [eV]

    % 1-D data interpolation
    for index = 1:length(Abs_array)
        Abs_array_eV{index} = interp1(E_A,fliplr(Abs_array{index}),E_S(1:end-1),'makima');
    end

    % create absorption table
    header = ["wavelength","energy"];
    data = [flipud(wavelength_S(1:end-1)),E_S(1:end-1)];
    for index = 1:length(Abs_array)
        header= horzcat(header, in.layers{index});
        data = horzcat(data, Abs_array_eV{index});
    end
    Abs_table = vertcat(header,data);
    save(append(in.prefix,"Abs.mat"), 'Abs_table');
    writematrix(Abs_table, append(in.prefix, "Abs.csv"), 'Delimiter',',');

    save(text);
end