function in=parameters()
%%reticolo input parameters
in.sym90 = false;       % true: Combination only when the x-axis and the y-axis have the same value

in.sym = true;         % true: The symmetry of the structure, more symmetry means shorter calculation time
in.pol = 0;             % polarization of the incident wave, TM pol=2  TE pol=0

in.npoints   = 9;       % point number of wavelength
in.lambdamin = 0.4;     % min wavelength
in.lambdamax = 1.2;     % max wavelength

%% Cross section along x0, y0, z0
in.x0 = 0;
in.y0 = 0;
in.z0 = 3.7;

%% x
% Number of Fourier terms
in.Mx = 0;
% period of nano structure
in.diam_x = 0.3;
% height of nano structure
in.height_nanostructure = 0.215;
% period of back grating
in.period_x = [2.4, 3.0];
% height of back grating
in.height_backgrating = 0.9;

in.My = in.Mx;                 % Number of Fourier terms in y
in.period_y = in.period_x;     % period in y
in.diam_y = in.diam_x;

