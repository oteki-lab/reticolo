function param_table=combination(in)
params = {};
keys = fieldnames(in);
for k =1:length(keys)
    key = char(keys(k));
    values = getfield(in, key);
    params = horzcat(params, values);
end
in_list = allcomb(params{:});
param_table = cell2table(table2cell(array2table(in_list)), 'VariableNames',keys);

