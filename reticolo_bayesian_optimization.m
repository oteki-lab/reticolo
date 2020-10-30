%% initialization
clear;
addpath(genpath('./'))
retio;
bayesian = py.importlib.import_module('bayesian');
py.importlib.reload(bayesian);

%% flags
notification    = false;    % true: send result mail (set address in sendMail.m)
cal_structure_x = false;    % true: calculate structure (x direction)
cal_structure_y = false;    % true: calculate structure (y direction)

take_over = '20201026135706568';     % directory name of simulation on the way
count_limit = 163;                   % limitation of iterative count, 0: only initialization

%% make output direcory in Results
dateString = datestr(datetime('now'),'yyyymmddHHMMssFFF');
disp(['make new directory: ',dateString]);
res_dir = ['Results\',dateString];
mkdir(res_dir);
mkdir([res_dir, '\graphs']);
mkdir([res_dir, '\bo_data']);
mkdir([res_dir, '\graphs/score']);
mkdir([res_dir, '\graphs/y_mean']);
mkdir([res_dir, '\graphs/acq']);

%% make inpout list & copy input file into output directory
mkdir([res_dir, '\input']);
copyfile('parameters.m', [res_dir,'\input\parameters.m'])
copyfile('structure.m', [res_dir,'\input\structure.m'])
copyfile('high_parameters.m', [res_dir,'\input\high_parameters.m'])
data = combination(parameters());
save([res_dir, '\input_list.mat'], 'data');

%% initizlize high parameters
hp_table_csv = [res_dir, '\hp_table.csv'];
hp_steps_csv = [res_dir, '\hp_steps.csv'];
hp_steps_mat = [res_dir, '\hp_steps.mat'];
[hp_table, hp_list] = high_parameters();
writetable(hp_table, hp_table_csv, 'Delimiter',',');

%% initilize steps
init_hp_steps_num = 0;
import_flag = false;
if isfile(['Results\', take_over, '\hp_steps.mat'])
    copyfile(['Results\', take_over, '\*'], [res_dir, '\']);
    load(hp_steps_mat)
    init_hp_steps_num = height(hp_steps);
    import_flag = true;
    count_limit = count_limit + init_hp_steps_num;
else
    hp_steps = [];
    params = {};
    %for i =1:height(hp_table)
    %    params = horzcat(params, [hp_table{i,2}, hp_table{i,3}]);
    %end
    params = horzcat(params, [0.10, 0.5, 1.0, 1.5, 2.5, 3.00]);
    params = horzcat(params, [0.05, 0.4, 0.9, 1.4, 2.4, 2.95]);
    init_steps = cell2table(table2cell(array2table(allcomb(params{:}))), 'VariableNames',hp_table{:,'name'});
    count_limit = count_limit + height(init_steps);
end

%% Run simulation
count = init_hp_steps_num;
%parpool
while(true)
    count = count + 1;
    disp(['Simulation: ', int2str(count), '/', int2str(count_limit)]);
    item_dir = append(res_dir, "\no_", int2str(count));
    mkdir(item_dir);
    attachments = strings(3:1);

    try
        % estimate next step by bayesian optimization
        if import_flag==false && count<=height(init_steps)
            next_step = table();
            for i=1:length(init_steps{count,:})
                key = char(init_steps.Properties.VariableNames(i));
                next_step.(key) = init_steps.(key)(count);
            end
        else
            next_step = struct2table(struct(py.bayesian.search(res_dir, hp_table_csv, hp_steps_csv)));
        end
        disp(next_step)
        % load input parameters
        in = table2struct(data(1,:));
        skip_flag = true;
        for i=1:length(next_step.Properties.VariableNames)
            key = char(next_step.Properties.VariableNames(i));
            in.(key) = next_step.(key);
        end

        % get structure parameters
        [in.params, in.layers] = structure(in);

        % output file name
        in.prefix = append(item_dir, "\no_", int2str(count), "_");
        in.res = ['period_',int2str(in.period_x*1000),'_diam_',int2str(in.diam_x*1000),'wav',int2str(in.lambdamin*1000),'_',int2str(in.lambdamax*1000),'_npoints',int2str(in.npoints),'_Fourier',int2str(in.Mx)];

        % calculate absorption
        %attachments(length(attachments)+1) = reticolo_eng(in);

        % calculate current density
        %J_table = current_density(in);
        %next_step.score = str2double(J_table(2,10))*1000;
        next_step.score = test(next_step.('Eg'), next_step.('Eic'));

        hp_steps = [hp_steps;next_step];
        writetable(hp_steps, hp_steps_csv, 'Delimiter',',');
        save(hp_steps_mat, 'hp_steps')
        disp(hp_steps)


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
    
    if count>=count_limit; break; end
end
delete(gcp('nocreate'))
clear;
