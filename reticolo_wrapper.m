addpath(genpath('./'))
clear;retio;

cal_absorption = true;      % ture: calculate absorption
cal_structure = true;       % ture: calculate structure

% make output direcory in Results
dateString = datestr(datetime('now'),'yyyyMMddHHmmssFFF');
disp(['make new directory: ',dateString]);
%res_dir = ['Results\',dateString];
res_dir = ['C:\Users\oteki\Dropbox\reticolo_data\',dateString];

mkdir(res_dir);

% copy input file
copyfile('input.py', [res_dir,'\input.py'])
copyfile('init_structure.m', [res_dir,'\init_structure.m'])

% make inpout list
command = ['python Codes\combination.py input.py ',res_dir];
status = dos(command);
data=load([res_dir,'\input_list']);
l = length(data.input_data);

% Run simulation
parpool
for index = 1:length(data.input_data)
    disp(['Simulation: ', int2str(index), '/', int2str(l)]);
    % load input parameters
    in = data.input_data{index};
    
    % get structure parameters
    [params, layers] = init_structure(in);
    
    % output file name
    %res = ['no_',int2str(index),'period_',int2str(in.period_x*1000),'_diam_',int2str(in.diam_x*1000),'wav',int2str(in.lambdamin*1000),'_',int2str(in.lambdamax*1000),'_npoints',int2str(in.npoints),'_Fourier',int2str(in.Mx)];
    res = ['period_',int2str(in.period_x*1000),'_diam_',int2str(in.diam_perc*1000),'_height_',int2str(in.height*1000),'.mat'];
    
    % calculate absorption
    if cal_absorption == true
        reticolo_eng(index, in, params, layers, res_dir, res);
    end
    
    % calculate structure
    if cal_structure == true
        in.Mx = 0;
        in.My = 0;
        in.npoints = 1;
        reticolo_eng(index, in, params, layers, res_dir, res);
    end
end
delete(gcp('nocreate'))
