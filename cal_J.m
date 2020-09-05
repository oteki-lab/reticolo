% constants
E_G=1.424;          % bandgap
h=6.63e-34;         % planck constant
c=3e8;              % speed of light
q=1.60e-19;         % elementary charge

% input file path
dir = 'results';
file = 'period_2400_diam_27wav400_1200_nbpoints101_Fourier0.mat';
text = [dir,'\',file];
number=1;

% read wavelength range & aborption rate
load(text,'Abs_array','layers','wavelength','lambdamin','lambdamax');
wavelength_A=wavelength;
wavmin=lambdamin*1000;
wavmax=lambdamax*1000;


load('Solarspectrum.mat') % spectra loaded in [W.m^-2.nm^-1]
% wavelength
wwavelength = wavelength(find(wavelength==wavmin):find(wavelength==wavmax));    % wavelength from wavmin to wavmax [nm]
E = flipud(h*c./(wwavelength*1e-9))/q;                                          % from wavelength [nm] to energy [eV]
% solar spectrum
AM1_5G = AM1_5G(find(wavelength==wavmin):find(wavelength==wavmax));             % solar spectrum from wavmin to wavmax [W.m^-2.nm^-1]
irr_nm = 1e-4*AM1_5G;                                                           % solar spectrum from [W.m^-2.nm^-1] to [W.cm^-2.nm^-1]
irr_eV = flipud(irr_nm*1e9)*h*c./(q*E.^2);                                      % solar spectrum from [W.cm-2.nm^-1] to [W.cm^-2.eV^-1]
J_tot = nansum((irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));           % Total current that can be converted from sunlight
% Heat irradiance
irr_kin_eV = max((1-E_G./E).*irr_eV,0);                                         % the kinetic power in [W.cm^-2.eV^-1]
irr_kin_nm = flipud(irr_kin_eV*q.*E.^2/(1e9*h*c));                              % the kinetic power from [W.cm-2.eV^-1] to [W.cm^-2.nm^-1]
P_kin_tot = nansum(irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));                 % Heat irradiance, or ideal thermalization intensity, assuming A=1, in W.cm^-2


% 1-D data interpolation
E_A=fliplr(h*c./(q*wavelength_A*1e-6));   % from wavelength [nm] to energy [eV]
E=fliplr(E);                            % Reverse order according to E_A
for index = 1:length(Abs_array)
    Abs_array_eV{index} = interp1(E_A,fliplr(Abs_array{index}),E,'makima');
end

% create absorption table
header = ["wavelength","energy"];
data = [flipud(wwavelength),E];
for index = 1:length(Abs_array)
    header= horzcat(header, layers{index});
    data = horzcat(data, Abs_array_eV{index});
end
Abs_table = vertcat(header,data);

% create J table
J_header = ["No"];
J_array = [number];
for index = 1:length(Abs_array)
    J_header = horzcat(J_header, layers{index});
    J_layer = nansum(Abs_array_eV{index}(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));
    J_array = horzcat(J_array, J_layer);
end
J_table = vertcat(J_header,J_array);

% calculate J of active region
Abs_layers = containers.Map(layers(:,1), Abs_array_eV);
A_active = Abs_layers(char('Active region'));
P_kin(number)=nansum(A_active(1:end-1).*irr_kin_eV(1:end-1).*(E(2:end)-E(1:end-1)));
J(number)=nansum(A_active(1:end-1).*(irr_eV(1:end-1)./E(1:end-1)).*(E(2:end)-E(1:end-1)));

save([dir,'\J_',file]);
