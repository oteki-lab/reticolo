addpath(genpath('./'))
clear;retio;

pole = 0;                   % polarization of the incident wave, TM pol=2  TE pol=0
cal_absorption = true;      % ture: calculate absorption
cal_structure_x = true;     % ture: calculate structure (x direction)
cal_structure_y = true;     % ture: calculate structure (y direction)

% make output direcory in Results
dateString = datestr(datetime('now'),'yyyyMMddHHmmssFFF');
disp(['make new directory: ',dateString]);
res_dir = ['Results\',dateString];
mkdir(res_dir);

% copy input file
copyfile('input.py', [res_dir,'\input.py'])
copyfile('init_structure.m', [res_dir,'\init_structure.m'])

% make inpout list
command = ['python Codes\combination.py input.py ',res_dir];
status = dos(command);
input_list=load([res_dir,'\input_list']);
l = length(input_list.data);

% Run simulation
parpool
for index = 1:length(input_list.data)
    disp(['Simulation: ', int2str(index), '/', int2str(l)]);

    % load input parameters
    in = input_list.data{index};
    
    % get structure parameters
    [params, layers] = init_structure(in);
    
    % output file name
    res = ['no_',int2str(index),'_period_',int2str(in.period_x*1000),'_diam_',int2str(in.diam_x*1000),'wav',int2str(in.lambdamin*1000),'_',int2str(in.lambdamax*1000),'_npoints',int2str(in.npoints),'_Fourier',int2str(in.Mx)];
    
    % calculate absorption
    if cal_absorption == true
        horizontal_label = 'x';
        reticolo_eng(index, in, pole, horizontal_label, params, layers, res_dir, res);
    end
    
    % calculate structure
    in.Mx = 0;
    in.My = 0;
    in.npoints = 1;
    % x direction
    if cal_structure_x == true
        horizontal_label = 'x';
        reticolo_eng(index, in, pole, horizontal_label, params, layers, res_dir, res);
    end
    % y direction
    if cal_structure_y == true
        horizontal_label = 'y';
        reticolo_eng(index, in, pole, horizontal_label, params, layers, res_dir, res);
    end
end
delete(gcp('nocreate'))
