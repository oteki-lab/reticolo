%% initialization
addpath(genpath('./'))
clear;retio;

%% make output direcory in Results
res_dir = mknewdir(['Results\', datestr(datetime('now'), 'yyyymmddHHMMSSFFF')]);
disp(['make new directory: ', res_dir]);

%% copy input file into output directory
filelist = {'parameters.m', 'layers.m', 'structure.m'};
for filename=filelist; copyfile(string(filename), append(res_dir, '\', string(filename))); end

%% make parameters list
params_list = combination(parameters());
save([res_dir, '\params_list.mat'], 'params_list');

%% Run simulation
parpool
for index = 1:height(params_list)
    disp(['Simulation: ', int2str(index), '/', int2str(height(params_list))]);
    attachments = strings(5:1);

    try
        % get input parameters
        in = table2struct(params_list(index,:));

        % get layers parameters
        [in.props, in.layers] = layers(in);

        % make output directory
        item_dir = mknewdir(append(res_dir, "\no_", int2str(index)));
        if in.out_I_map; in.map_dir = mknewdir(append(item_dir,'\I_map')); end
        if in.trace_champ; in.cs_dir = mknewdir(append(item_dir,'\cross_section')); end

        % set output file name
        in.prefix = append(item_dir, "\no_", int2str(index), "_");
        in.res = ['period_',int2str(in.period_x*1000),'_diam_',int2str(in.diam_x*1000),'wav',int2str(in.lambdamin*1000),'_',int2str(in.lambdamax*1000),'_npoints',int2str(in.npoints),'_Fourier',int2str(in.Mx)];

        % calculate absorption
        if in.cal_absorption
            attachments(length(attachments)+1) = reticolo_eng(in);

            % calculate current density
            if in.cal_current
                current_density(in);
            end
        end

        % calculate structure
        if in.cal_structure
            in.cal_absorption = false;
            in.cal_field      = false;
            in.trace_champ    = true;
            in.Mx = 0;
            in.My = 0;
            in.npoints = 1;
            in.cs_dir = mknewdir(append(item_dir,'\structure'));
            reticolo_eng(in);
        end
        
       %%
        msg = "reticolo Simulation Done.";
        
    catch e
        disp(e)
        msg = "reticolo Simulation Error.";
        delete(gcp('nocreate'))
    end
    if in.notification; sendMail(msg, attachments); end
end
delete(gcp('nocreate'))
