addpath(genpath('./'))
clear;retio;
cal_absorption = true;
cal_structure = true;

% make output direcory in Results
dateString = datestr(datetime('now'),'yyyyMMddHHmmssFFF');
disp(['make new directory: ',dateString]);
res_dir = ['Results\',dateString];
mkdir(res_dir);

% make input file
copyfile('input.py', [res_dir,'\input.py'])
copyfile('input.py', [res_dir,'\init_structure.m'])
command = ['python Codes\combination.py input.py ',res_dir];
status = dos(command);

% Run simulation
data=load([res_dir,'\input_list']);
l = length(data.input_data);

parpool
for index = 1:length(data.input_data)
    disp(['Simulation: ', int2str(index), '/', int2str(l)]);
    in = data.input_data{index};
    
    res = ['no_',int2str(index),'period_',int2str(in.period_x*1000),'_diam_',int2str(in.diam_x*1000),'wav',int2str(in.lambdamin*1000),'_',int2str(in.lambdamax*1000),'_npoints',int2str(in.npoints),'_Fourier',int2str(in.Mx)];
    [params, layers] = init_structure(in);
    
    if cal_absorption == true
        reticolo_eng(index, in, params, layers, res_dir, res);
    end
    
    if cal_structure == true
        in.npoints = 1;
        reticolo_eng(index, in, params, layers, res_dir, res);
    end
end
delete(gcp('nocreate'))
