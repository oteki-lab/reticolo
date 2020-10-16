function in=parameters()
%%reticolo input parameters
in.asymmetry = true;    % True: Combination only when the x-axis and the y-axis have the same value

in.pol = 0;             % polarization of the incident wave, TM pol=2  TE pol=0
in.npoints = 16;        % point number of wavelength
in.lambdamin = 0.9;     % min wavelength
in.lambdamax = 1.2;     % max wavelength

%% x
% Number of Fourier terms
in.Mx = 0;
% period of nano structure
in.diam_x = 0.215;
% height of nano structure
in.height_nanostructure = 0.4;
% period of back grating
in.period_x = 2.4;
% height of back grating
in.height_backgrating = 0.9;

in.My = in.Mx;                 % Number of Fourier terms in y
in.period_y = in.period_x;     % period in y
in.diam_y = in.diam_x;
