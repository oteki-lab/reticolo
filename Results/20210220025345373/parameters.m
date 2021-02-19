function in=parameters()
%% flags
in.notification     = false;    % true: send result mail (set address in sendMail.m)
in.cal_absorption   = true;    % true: calculate absorption
in.cal_field        = true;    % true: calculate field intensity (3D)
in.out_I_map        = true;    % true: output field intensity (3D) data
in.cal_current      = true;    % true: calculate current density from absorption
in.trace_champ      = true;    % true: calculates a cross-section of the field intensity in each wavelength
in.cal_structure    = true;    % true: calculate a cross-section of the structure
in.cal_structure_yz = true;    % true: calculate a cross-section of the field intensity & structure (x direction)
in.cal_structure_xz = false;    % true: calculate a cross-section of the field intensity & structure (y direction)
in.cal_structure_xy = false;    % true: calculate a cross-section of the field intensity & structure (z direction)

in.sym90 = true;       % true: Combination only when the x-axis and the y-axis have the same value
in.sym = true;         % true: The symmetry of the structure, more symmetry means shorter calculation time

%% reticolo input parameters
in.pol       = 0;      % polarization of the incident wave, TM pol=2  TE pol=0
in.npoints   = 31;     % point number of wavelength
in.lambdamin = 0.90;   % min wavelength
in.lambdamax = 1.20;   % max wavelength

%% Cross section along x0, y0, z0
in.x0 = 0;
in.y0 = 0;
in.z0 = 1.310;

%% property parameters
% Number of Fourier terms
in.Mx = 0;
% height of nano structure
in.height_nanostructure = 0.215;
% height of back grating
in.height_backgrating = 0.9;
% period of back grating
in.period_x = [2.4];
% period of nano structure
in.diam_x = 0.0;    %0.3;

in.My = in.Mx;                 % Number of Fourier terms in y
in.period_y = in.period_x;     % period in y
in.diam_y = in.diam_x;

