%% initialization
addpath(genpath('./'))
clear;retio;

%% flags
notification    = false;    % true: send result mail (set address in sendMail.m)
cal_absorption  = true;     % true: calculate absorption
cal_structure_x = true;     % true: calculate structure (x direction)
cal_structure_y = true;     % true: calculate structure (y direction)

%% make output direcory in Results
dateString = datestr(datetime('now'),'yyyyMMddHHmmssFFF');
disp(['make new directory: ',dateString]);
res_dir = ['Results\',dateString];
mkdir(res_dir);

%% copy input file into output directory
copyfile('parameters.m', [res_dir,'\parameters.m'])
copyfile('structure.m', [res_dir,'\structure.m'])

%% make inpout list
data = combination(parameters());
save([res_dir, '\input_list.mat'], 'data');
l = height(data);

%% Run simulation
%parpool
for index = 1:l
    disp(['Simulation: ', int2str(index), '/', int2str(l)]);
    attachments = strings(3:1);

    try
        % load input parameters
        in = table2struct(data(index,:));

        % get structure parameters
        [in.params, in.layers] = structure(in);

        % output file name
        in.prefix = append(res_dir, "\no_", int2str(index), "_");
        in.res = ['period_',int2str(in.period_x*1000),'_diam_',int2str(in.diam_x*1000),'wav',int2str(in.lambdamin*1000),'_',int2str(in.lambdamax*1000),'_npoints',int2str(in.npoints),'_Fourier',int2str(in.Mx)];

        % calculate absorption
        if cal_absorption
            attachments(length(attachments)+1) = reticolo_eng(in);
        end

        % calculate structure
        in.Mx = 0;
        in.My = 0;
        in.npoints = 1;
        % x direction
        if cal_structure_x
            in.horizontal_label = 'x';
            attachments(length(attachments)+1) = reticolo_eng(in);
        end
        % y direction
        if cal_structure_y
            in.horizontal_label = 'y';
            attachments(length(attachments)+1) = reticolo_eng(in);
        end
        
        msg = "reticolo Simulation Done.";
        
    catch e
        disp(e)
        msg = "reticolo Simulation Error.";
        delete(gcp('nocreate'))
    end
    if notification; sendMail(msg, attachments); end
end
delete(gcp('nocreate'))
