

%%%%                Calculate diffraction in a system composed of 1-12 layers
%%%%                                         (17/06/2014)
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

clear;retio;
[prv,vmax]=retio([],inf*1i);

%% Parameters of the structure and the calculation


%%%%%% Wavelengths and angle of incidence
npoints=101;  % 1 for only structure
lambdamin=0.4;
lambdamax=1.2;
wavelength=linspace(lambdamin,lambdamax,npoints);
theta=[0,0];                          %angle of incidence in degrees

%%%%%% Geometric parameters
periodicity_x=2.4;                  % period in x
periodicity_y=periodicity_x;        % period in y
Nb_couches=18; %Number of layers
diam=0.215/4;

%diameter of each layer
Diameter_x=[
    diam*1,
    diam*2,
    diam*3,
    periodicity_x,
    periodicity_x,
    periodicity_x,
    periodicity_x,
    periodicity_x*10/10,
    periodicity_x*9/10,
    periodicity_x*8/10,
    periodicity_x*7/10,
    periodicity_x*6/10,
    periodicity_x*5/10,
    periodicity_x*4/10,
    periodicity_x*3/10,
    periodicity_x*2/10,
    periodicity_x*1/10,
    diam
];

% Thicknesses, from top to bottom   (0 si if no layer)
Height=[
    0.4/3,                          % AlInP
    0.4/3,                          % AlInP
	0.4/3,                          % AlInP
	0.04,                           % AlInP
	0.16,                           % GaAs
	0.14,                           % QD
	1.7,                            % GaAs
	0.04,                           % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.9/9,                          % AlInP
	0.01                            % Ag
];

%%%%%% Refraction indices (from top to bottom), can be a function of the wavelength
nh=1;

n = [
    ones(size(wavelength)),                  % AlInP
	ones(size(wavelength)),                  % AlInP
    ones(size(wavelength)),                  % AlInP
    ones(size(wavelength)),                  % AlInP
    retindice_chen(wavelength,4.707),        % emitter GaAs
    retindice_chen(wavelength,4.708),        % QD Al80Ga20As 
    retindice_chen(wavelength,4.707),        % base GaAs
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    1.58*ones(size(wavelength)),             % AlInP
    retindice_chen(wavelength,1.72);         % Ag
];
nm = [
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    0.00*ones(size(wavelength)),            %
    0.00*ones(size(wavelength)),            %
    0.00*ones(size(wavelength)),            %
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    retindice_chen(wavelength,4.802),       % 
    0.00*ones(size(wavelength))             %
];

nsub=retindice_chen(wavelength,1.72);       % Ag argent index of the substrate

%%%%%% Numerical parameters
pol=0;                              % polarization of the incident wave, TM pol=2  TE pol=0
                                    % For normal incidence, TM <=> H//y and TE <=> E//y
sym=[pol-1,pol-1,0,0];              % The symmetry of the structure, more symmetry means shorter calculation time
% IMPORTANT: To be changed if non-normal incident or if non-rectangular
% structures
% if theta(1)==0 && theta(2)==0;sym=[pol-1,pol-1,0,0];end;
% if theta(1)~=0 && theta(2)==0;sym=[0,pol-1,0,0];end;
% if theta(1)==0 && theta(2)~=0;sym=[1-pol,0,0,0];end;
% if theta(1)~=0 && theta(2)~=0;sym=[];end;

%% 
Mx=0;                                % Number of Fourier terms in x
My=Mx;                               % Number of Fourier terms in y
op_granet=0;                         % If 1, RCWA is modified to improve convergence (Transforms the real coordinates at discontinuities)
% IMPORTANT: this parameter is tricky to use, and does not work out of normal incidence. Better keep it at zero

cal_abs=1;                           % If 1, calculate absorption in each layer
Nb_pts_z=10;                         % Number of points in z to calculate absorption, when absorption is calculated in each layer

cal_champ=0;                         % If 1, calculate the field in layer N_semicon
% IMPORTANT: only one wavelength (calculate can be quite heavy, depending on Nb_pts_z_semicon)
N_semicon=3;                         % Layer where field is calculated (if cal_champ=1)
Nb_pts_z_semicon=50;                 % Number of points in the z direction to calculate the field in N_semicon (if cal_champ=1)

trace_champ=0;                       % si 1,calculates a cross-section of the field
% IMPORTANT: only one wavelength
x0=0;                                % Cross section along x=x0 if trace_champ=1 ([] if the cross-section is along another direction)
y0=[];                               % Cross section along y=y0 if trace_champ=1 ([] if the cross-section is along another direction)
z0=[];                               % Cross section along z=z0 if trace_champ=1 ([] if the cross-section is along another direction)
% note: z=0 corresponds to the bottom of the considered stack, at a depth h_sub inside the substrate
h_air=0.05;                          % Thickness in incident medium to represent the cross-section (trace_champ=1)
h_sub=0.05;                          % Thickness in the substrate to represent the cross-section  (trace_champ=1)
h_2pts=20;                           % distance between 2 points in z pour the layers whose thickness is higher than 1??m (trace_champ=1)
op_objet=0;                          % If 1, plot the geometry to verify the calculated structure is correct

if op_granet==1
    Bx=500;Ax=0.02/Bx;By=Bx;Ay=Ax;
    diameter_x=Diameter_x(1:(Nb_couches));diameter_y=diameter_x;
    xdisc=[-diameter_x/2,diameter_x/2];ydisc=[-diameter_y/2,diameter_y/2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code (for parallel computing of the absorption)

% Initialisation of the parameters for parallel computing
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
H=Height(1:(Nb_couches));
if cal_abs==1||cal_champ==1||trace_champ==1;op_retcouche=1; else; op_retcouche=0; end
if H(Nb_couches)<1e-5;disp('WARNING : There is a problem in the definition of the layers number !!'); return; end
if trace_champ==1&&isempty(x0)==1&&isempty(y0)==1&&isempty(z0)==1;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end
if trace_champ==1&&isempty(x0)==0&&isempty(y0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end
if trace_champ==1&&isempty(x0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end
if trace_champ==1&&isempty(y0)==0&&isempty(z0)==0;disp('WARNING : There is a problem in the definition of the desired cross section for plotting the field (trace_champ=1) !!'); return; end

parpool
parfor zou=1:length(wavelength)
%for zou=1:length(wavelength)
    disp(['Calculation n-' int2str(zou) ' of ' int2str(length(wavelength))])
    
    inc=[];
    sym=[];
    e0=[];
    o0=[];
    % R0_TE_TE_vect=zeros(1);
    % R0_TE_TM_vect=R0_TE_TE_vect;
    % R0_TM_TM_vect=R0_TE_TE_vect;
    % R0_TM_TE_vect=R0_TE_TE_vect;
    ref_TE_TE_vect=zeros(1);ref_TE_TM_vect=ref_TE_TE_vect;ref_TM_TM_vect=ref_TE_TE_vect;ref_TM_TE_vect=ref_TE_TE_vect;
    A_tot_vect=zeros(1);
    A_sub_vect=zeros(1);
    Abs_vect=zeros(Nb_couches,1);
    Abs_plots_vect=zeros(Nb_couches,1);
    test_vect=zeros(1);
    Einc=[];Hinc=[];E_semicon=[];H_semicon=[];x_semicon=[];y_semicon=[];z_semicon=[];W_semicon=[];
    Ex=[];Ey=[];Ez=[];Hx=[];Hy=[];Hz=[];
    
    xx=[];yy=[];zz=[];indice=[];
    lwa=wavelength(zou);
    k0=2*pi/lwa;
    %tic
    period=[periodicity_x,periodicity_y];

    ns=nsub(zou);
    Number = [];
    for index = 1:length(n(:,1))
        if length(n(index,:))==1; Number=[Number,n(index)]; else; Number=[Number,n(index, zou)]; end
    end
    
    Numberm = [];
    for index = 1:length(nm(:,1))
        if length(nm(index,:))==1; Numberm=[Numberm,nm(index)]; else; Numberm=[Numberm,nm(index, zou)]; end
    end
    
    N=Number(1:Nb_couches);
    Nm=Numberm(1:Nb_couches);
    diameter_x=Diameter_x(1:(Nb_couches));
    diameter_y=diameter_x;
    kx=k0*nh*sin(theta(1)*pi/180);
    ky=k0*nh*sin(theta(2)*pi/180);
    beta=[kx,ky];
    
    uh=retu(period,{nh,k0});
    ub=retu(period,{ns,k0});
    
    if op_granet==1
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym,{[],{xdisc,Ax,Bx,ydisc,Ay,By}});ah=retcouche(init,uh,1);
    else
        init=retinit(period,[-Mx,Mx,-My,My],beta,sym);ah=retcouche(init,uh,op_retcouche);
    end
    ab=retcouche(init,ub,op_retcouche);
    
    u=[];
    a=[];
    for az=1:Nb_couches
        if Nm(az)==0
            u{az}=retu(period,{N(az),k0});
        else
            % u{az}=retu(period,{N(az),[0,0,diameter_x,diameter_y,Nm(az),Ntre],[-diameter_x/2+w_rectangle/2,diameter_y/2+h_rectangle/2,w_rectangle,h_rectangle,Nm(az),Ntre],[diameter_x/2+h_rectangle/2,-diameter_y/2+w_rectangle/2,h_rectangle,w_rectangle,Nm(az),Ntre ],k0});
            if az==1||az==2||az==3
                structure1=[-1200*10/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure2=[-1200*10/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure3=[-1200*10/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure4=[-1200*10/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure5=[-1200*10/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure6=[-1200*10/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure7=[-1200*10/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure8=[-1200*10/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure9=[-1200*10/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure10=[-1200*10/11,1200*8,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure11=[-1200*10/11,1200*10,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure12=[-1200*8/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure13=[-1200*8/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure14=[-1200*8/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure15=[-1200*8/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure16=[-1200*8/10,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure17=[-1200*8/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure18=[-1200*8/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure19=[-1200*8/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure20=[-1200*8/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure21=[-1200*8/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure22=[-1200*8/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure23=[-1200*6/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure24=[-1200*6/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure25=[-1200*6/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure26=[-1200*6/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure27=[-1200*6/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure28=[-1200*6/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure29=[-1200*6/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure30=[-1200*6/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure31=[-1200*6/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure32=[-1200*6/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure33=[-1200*6/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure34=[-1200*4/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure35=[-1200*4/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure36=[-1200*4/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure37=[-1200*4/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure38=[-1200*4/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure39=[-1200*4/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure40=[-1200*4/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure41=[-1200*4/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure42=[-1200*4/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure43=[-1200*4/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure44=[-1200*4/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure45=[-1200*2/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure46=[-1200*2/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure47=[-1200*2/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure48=[-1200*2/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure49=[-1200*2/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure50=[-1200*2/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure51=[-1200*2/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure52=[-1200*2/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure53=[-1200*2/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure54=[-1200*2/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure55=[-1200*2/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure56=[0,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure57=[0,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure58=[0,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure59=[0,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure60=[0,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure61=[0,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure62=[0,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure63=[0,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure64=[0,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure65=[0,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure66=[0,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure67=[1200*2/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure68=[1200*2/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure69=[1200*2/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure70=[1200*2/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure71=[1200*2/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure72=[1200*2/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure73=[1200*2/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure74=[1200*2/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure75=[1200*2/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure76=[1200*2/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure77=[1200*2/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure78=[1200*4/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure79=[1200*4/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure80=[1200*4/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure81=[1200*4/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure82=[1200*4/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure83=[1200*4/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure84=[1200*4/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure85=[1200*4/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure86=[1200*4/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure87=[1200*4/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure88=[1200*4/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure89=[1200*6/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure90=[1200*6/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure91=[1200*6/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure92=[1200*6/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure93=[1200*6/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure94=[1200*6/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure95=[1200*6/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure96=[1200*6/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure97=[1200*6/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure98=[1200*6/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure99=[1200*6/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure100=[1200*8/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure101=[1200*8/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure102=[1200*8/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure103=[1200*8/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure104=[1200*8/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure105=[1200*8/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure106=[1200*8/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure107=[1200*8/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure108=[1200*8/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure109=[1200*8/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure110=[1200*8/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                structure111=[1200*10/11,-1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure112=[1200*10/11,-1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure113=[1200*10/11,-1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure114=[1200*10/11,-1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure115=[1200*10/11,-1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure116=[1200*10/11,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure117=[1200*10/11,1200*2/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure118=[1200*10/11,1200*4/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure119=[1200*10/11,1200*6/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure120=[1200*10/11,1200*8/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                structure121=[1200*10/11,1200*10/11,diameter_x(az),diameter_y(az),Nm(az),Ntre];

                
                texture={N(az),structure1,structure2,structure3,structure4,structure5,structure6,structure7,structure8,structure9,structure10,structure11,structure12,structure13,structure14,structure15,structure16,structure17,structure18,structure19,structure20,structure21,structure22,structure23,structure24,structure25,structure26,structure27,structure28,structure29,structure30,structure31,structure32,structure33,structure34,structure35,structure36,structure37,structure38,structure39,structure40,structure41,structure42,structure43,structure44,structure45,structure46,structure47,structure48,structure49,structure50,structure51,structure52,structure53,structure54,structure55,structure56,structure57,structure58,structure59,structure60,structure61,structure62,structure63,structure64,structure65,structure66,structure67,structure68,structure69,structure70,structure71,structure72,structure73,structure74,structure75,structure76,structure77,structure78,structure79,structure80,structure81,structure82,structure83,structure84,structure85,structure86,structure87,structure88,structure89,structure90,structure91,structure92,structure93,structure94,structure95,structure96,structure97,structure98,structure99,structure100,structure101,structure102,structure103,structure104,structure105,structure106,structure107,structure108,structure109,structure110,structure111,structure112,structure113,structure114,structure115,structure116,structure117,structure118,structure119,structure120,structure121,k0};
                
                
            else
                structure=[0,0,diameter_x(az),diameter_y(az),Nm(az),Ntre];
                texture={N(az),structure,k0};
                
                
            end
            
            u{az}=retu(period,texture);
        end
        a{az}=retcouche(init,u{az},op_retcouche);
    end
    struct_test=[];
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
        [X,wX]=retgauss(-periodicity_x/2,periodicity_x/2,15,10,Xdisc);
        [Y,wY]=retgauss(-periodicity_y/2,periodicity_y/2,15,10,Ydisc);
        [x,y]=retgranet(init,'num2phys',X,Y);
        [X,Y,Xder,Yder]=retgranet(init,x,y);
        wx=wX./Xder;wy=wY./Yder;
    else
        [x,wx]=retgauss(-periodicity_x/2,periodicity_x/2,15,10,[-unique(diameter_x)/2,unique(diameter_x)/2]);
        [y,wy]=retgauss(-periodicity_y/2,periodicity_y/2,15,10,[-unique(diameter_y)/2,unique(diameter_y)/2]);
    end
    
    if cal_abs==1||cal_champ==1||trace_champ==1
        sb_norm=retb(init,ah,-0.1,0,[],[]);
        tab_norm=[0,1,1];
        if size(ef.inc.teta,2)==1
            inc=1;
        elseif size(ef.inc.teta,2)==2
            inc=[0,0];
            if pol==0;inc(ef.inc.TE)=1; elseif pol==2; inc(ef.inc.TM)=1; end
        end
        for az=1:(Nb_couches)
            [ei,zi]=retchamp(init,{ah},sh,sb_norm,inc,{x,y},tab_norm,[],1:6,1,1,1:6);
            flux=retpoynting(ei,[0,0,-1],wx,wy,[]);
            inc=1/sqrt(flux).*inc;

            [einc,zinc]=retchamp(init,{ah},sh,sb_norm,inc,{0,0},tab_norm,[],1:6,1,1,1:6);
            Einc=squeeze(einc(:,:,:,1:3));
            Hinc=squeeze(einc(:,:,:,4:6));
        end
    end
    
    if cal_abs==1
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
        
        tab=[];tab2=tab;
        struct={ab};
        for az=1:Nb_couches
            tab=[tab;[H(az),az+1,0]];
            tab2=[tab2;[0,az+1,1];[H(az),az+1,0]];
            struct={struct{:},a{az}};
        end
        tab(num,3)=Nb_pts_z;
        tab2=[tab2;[0,1,1];[1,1,0];[0,1,1]];
        
        [e,z,wz,o]=retchamp(init,struct,sh,sb,inc,{x,y},tab,[],[1:6]+7.25i,1,1,1:6);
        [e2,z2,wz2,o2]=retchamp(init,struct,sh,sb,inc,{x,y},tab2,[],[1:6]+7.25i,1,1,1:6);
        for ii=1:3
            o(:,:,:,ii+3)=o(:,:,:,ii+3)./o(:,:,:,ii);
            o(:,:,:,ii)=1;
        end
        [Wz,Wx,Wy]=ndgrid(wz,wx,wy);
        W=Wz.*Wx.*Wy;
        flux_poyn=retpoynting(e2,[0,0,-1],wx,wy,[]);
        
        A_sub_vect=flux_poyn(2)-flux_poyn(1);
        
        for az=1:Nb_couches
            Abs_vect(Nb_couches-az+1)=flux_poyn(az+2)-flux_poyn(az+1);
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
        test(:,zou)=test_vect;
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
        
        [e_semicon,z_semicon,wzs,o_semicon]=retchamp(init,struct,sh,sb,inc,{x,y},tab_semicon,[],[1:6]+7.25i,1,1,1:6);
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
            [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,y0},tab0,[],[1:6]+7.25i,1,1,1:6);
        elseif isempty(y0)==1&&isempty(z0)==1
            [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{x0,yy},tab0,[],[1:6]+7.25i,1,1,1:6);
        elseif isempty(x0)==1&&isempty(y0)==1
            tab0(:,3)=0;
            HH=cumsum(tab0(:,1));
            numz=1;
            while (HH(end)-z0)>HH(numz);numz=numz+1; end
            tab1=[tab0(1:numz-1,:);[tab0(numz,1)-(z0-sum(tab0(numz+1:end,1))),numz,0];[0,numz,1];[z0-sum(tab0(numz+1:end,1)),numz,0];tab0(numz+1:end,:)];
            [e0,zz,wz,o0]=retchamp(init,struct0,sh,sb,inc,{xx,yy},tab1,[],[1:6]+7.25i,1,1,1:6);
        end
        
        for ii=1:3
            o0(:,:,:,ii+3)=o0(:,:,:,ii+3)./o0(:,:,:,ii);
            o0(:,:,:,ii)=1;
        end
        indice=squeeze(sqrt(o0(:,:,:,4)));
        Ex=squeeze(e0(:,:,:,1));Ey=squeeze(e0(:,:,:,2));Ez=squeeze(e0(:,:,:,3));
        Hx=squeeze(e0(:,:,:,4));Hy=squeeze(e0(:,:,:,5));Hz=squeeze(e0(:,:,:,6));
    end

end
delete(gcp('nocreate'))
%% Saving and plotting output data

%%%% Example to plot a cross section with trace_champ=1, x0=[], y0=0, z0=[]
%%%% Hy as a function of x and z
%%%% The separation between the layers is in white (indice contains the position of all indices)
if trace_champ == 1
    figure
    pcolor(xx,zz,abs(Ey).^2);
    shading interp;
    colormap(jet);
    hold on
    contour(real(xx),real(zz),real(indice),'black','linewidth',2);
    xlabel('x')
    ylabel('z')
end

if cal_abs == 1
    %%%% Example to save the data into a file
    text=['period_',int2str(periodicity_x*1000),'_diam_',int2str(diam*1000),'wav',int2str(wavelength(1)*1000),'_',int2str(wavelength(end)*1000),'_nbpoints',int2str(length(wavelength)),'_Fourier',int2str(Mx),'.mat'];
    save(text);

    %%%% Example to plot the absorption
    figure
    %plot(wavelength,A_tot(1,:),'Linewidth',3)
    %legend({'Total'})
    plot(wavelength,Abs(1,:)+Abs(2,:)+Abs(3,:)+Abs(4,:), wavelength,Abs(5,:), wavelength,Abs(6,:), wavelength,Abs(7,:), wavelength,Abs(8,:)+Abs(9,:)+Abs(10,:)+Abs(11,:)+Abs(12,:)+Abs(13,:)+Abs(14,:)+Abs(15,:)+Abs(16,:)+Abs(17,:), wavelength,Abs(18,:), wavelength,Abs(5,:)+Abs(6,:)+Abs(7,:), 'Linewidth',3)
    legend({'AlInP(1-4)','GaAs(5)','QD(6)','GaAs(7)','AlInP(8-17)','Silver mirror(18)','active region(4-6)'})
    xlabel('\lambda (μm)')
    ylabel('Absorption')
    xlim([min(wavelength) max(wavelength)])
    ylim([0 1])
    set(gca,'Fontsize',12)
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gcf,'color','w');
    box on
end