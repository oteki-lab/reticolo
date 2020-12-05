function y=current_density(in)
%% Saving and plotting output data
%%%% Save the data into a file
text=append(in.prefix, in.res, '.mat');
if exist(text, 'file')
    load(text);
    
    J_sum_table = ["No"];
    for index = 1:length(Abs_array)
        J_sum_table = horzcat(J_sum_table, in.layers{index});
    end
 
    J_table = ["wavelength (nm)"];
    for index = 1:length(Abs_array)
        J_table = horzcat(J_table, in.layers{index});
    end
    
    % solar spectrum
    AM_list = ["AM1.5G","AM1.5D","AM0"];
    load('Solarspectrum.mat','AM1_5G','AM1_5D','AM0')                                           % spectra loaded in [W.m^-2.nm^-1]
    AM1_5G = AM1_5G(w_range);                                                                   % solar spectrum from wavmin to wavmax [W.m^-2.nm^-1]
    AM1_5D = AM1_5D(w_range);                                                                   % solar spectrum from wavmin to wavmax [W.m^-2.nm^-1]
    AM0 = AM0(w_range);                                                                         % solar spectrum from wavmin to wavmax [W.m^-2.nm^-1]
    AM = {AM1_5G, AM1_5D, AM0};
    
    for i=1:length(AM)
        irr_nm = 1e-4*AM{i};                                                                    % solar spectrum from [W.m^-2.nm^-1] to [W.cm^-2.nm^-1]
        irr_eV = flipud(irr_nm*1e9)*h*c./(q*E_S.^2);                                            % solar spectrum from [W.cm-2.nm^-1] to [W.cm^-2.eV^-1]
        J_tot = 1e3*sum((irr_eV(1:end-1)./E_S(1:end-1)).*(E_S(2:end)-E_S(1:end-1)),'omitnan');  % Total current that can be converted from sunlight
        % Heat irradiance
        irr_kin_eV = max((1-E_G./E_S).*irr_eV,0);                                               % the kinetic power in [W.cm^-2.eV^-1]
        irr_kin_nm = flipud(irr_kin_eV*q.*E_S.^2/(1e9*h*c));                                    % the kinetic power from [W.cm-2.eV^-1] to [W.cm^-2.nm^-1]
        P_kin_tot = 1e3*sum(irr_kin_eV(1:end-1).*(E_S(2:end)-E_S(1:end-1)),'omitnan');          % Heat irradiance, or ideal thermalization intensity, assuming A=1, in W.cm^-2

        % create J table
        J_sum_array = [AM_list(i)];
        J_array = [wavelength_S(1:end-1)];
        for index = 1:length(Abs_array)
            J = 1e3*Abs_array_eV{index}(1:end).*(irr_eV(1:end-1)./E_S(1:end-1)).*(E_S(2:end)-E_S(1:end-1));
            J_layer = sum(J,'omitnan');
            J_sum_array = horzcat(J_sum_array, J_layer);
            J_array = horzcat(J_array, flip(J));
        end
        J_sum_table = vertcat(J_sum_table, J_sum_array);
        J_table_tmp = vertcat(J_table, J_array); 
        filename = append(in.prefix, "J_", AM_list(i), ".csv");
        writematrix(J_table_tmp, filename, 'Delimiter',',');
        J_array=[];
    end
    filename = append(in.prefix,"J.mat");
    save(filename, 'J_sum_table');
    filename = append(in.prefix, "J.csv");
    writematrix(J_sum_table, filename, 'Delimiter',',');

    save(text);
else
    disp('Not found absorption data');
end
