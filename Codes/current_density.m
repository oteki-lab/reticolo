function y=current_density(in)
%% Saving and plotting output data
%%%% Save the data into a file
text=append(in.prefix, in.res, '.mat');
if exist(text, 'file')
    load(text);
    
    % solar spectrum
    load('Solarspectrum.mat','AM1_5G','AM1_5D','AM0')                               % spectra loaded in [W.m^-2.nm^-1]
    AM1_5G = AM1_5G(w_range);                                                       % solar spectrum from wavmin to wavmax [W.m^-2.nm^-1]
    irr_nm = 1e-4*AM1_5G;                                                           % solar spectrum from [W.m^-2.nm^-1] to [W.cm^-2.nm^-1]
    irr_eV = flipud(irr_nm*1e9)*h*c./(q*E_S.^2);                                    % solar spectrum from [W.cm-2.nm^-1] to [W.cm^-2.eV^-1]
    J_tot = nansum((irr_eV(1:end-1)./E_S(1:end-1)).*(E_S(2:end)-E_S(1:end-1)));     % Total current that can be converted from sunlight
    % Heat irradiance
    irr_kin_eV = max((1-E_G./E_S).*irr_eV,0);                                       % the kinetic power in [W.cm^-2.eV^-1]
    irr_kin_nm = flipud(irr_kin_eV*q.*E_S.^2/(1e9*h*c));                            % the kinetic power from [W.cm-2.eV^-1] to [W.cm^-2.nm^-1]
    P_kin_tot = nansum(irr_kin_eV(1:end-1).*(E_S(2:end)-E_S(1:end-1)));             % Heat irradiance, or ideal thermalization intensity, assuming A=1, in W.cm^-2

    % create J table
    J_header = ["No"];
    J_array = ["AM1.5"];
    for index = 1:length(Abs_array)
        J_header = horzcat(J_header, in.layers{index});
        J_layer = nansum(Abs_array_eV{index}(1:end-1).*(irr_eV(1:end-1)./E_S(1:end-1)).*(E_S(2:end)-E_S(1:end-1)));
        J_array = horzcat(J_array, J_layer);
    end
    J_table = vertcat(J_header,J_array);
    filename = append(in.prefix,"J.mat");
    save(filename, 'J_table');

    % calculate J of active region
    Abs_layers = containers.Map(in.layers(:,1), Abs_array_eV);
    A_active = Abs_layers(char('Active region'));
    P_kin=nansum(A_active(1:end-1).*irr_kin_eV(1:end-1).*(E_S(2:end)-E_S(1:end-1)));
    J=nansum(A_active(1:end-1).*(irr_eV(1:end-1)./E_S(1:end-1)).*(E_S(2:end)-E_S(1:end-1)));
    save(text);
else
    disp('Not found absorption data');
end
