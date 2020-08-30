addpath(genpath('./'))
clear;retio;

command = "python main.py";
status = dos(command);

data=load('input');
l = length(data.input_data);
count = 0;
for index = 1:length(data.input_data)
    count = count + 1;
    disp(['Simulation: ', int2str(count), '/', int2str(l)]);
    in = data.input_data{index};
    [params, layers, res] = init_structure(in);
    reticolo_eng(in, params, layers, res);
end
