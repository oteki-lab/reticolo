function [hp_table, hp_list]=high_parameters()
%high_parameters = {
%    {'period_x', 1.0, 3.0, 0.1}, {'period_y', 1.0, 3.0, 0.1}
%};

%% set high parameters (name, min, max, step)
high_parameters = {
    %'period_x' 1.0 3.0 0.1;
    'h2' 0.0 10.0 0.5;
    'h' 0.0 10.0 0.5
};

%% shape high parameters as table
hp_table = cell2table(table2cell(array2table(high_parameters)), 'VariableNames',{'name','min','max','step'});

%% list up combination of all high parameters
params = {};
for i =1:height(hp_table)
    params = horzcat(params, (hp_table{i,2}:hp_table{i,4}:hp_table{i,3}));
end
hp_list = cell2table(table2cell(array2table(allcomb(params{:}))), 'VariableNames',hp_table{:,'name'});
