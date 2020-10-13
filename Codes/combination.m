function param_table=combination(in)
params = {};
keys = fieldnames(in);
for k =1:length(keys)
    key = char(keys(k));
    values = in.(key);
    params = horzcat(params, values);
end
in_list = allcomb(params{:});
t = cell2table(table2cell(array2table(in_list)), 'VariableNames',keys);

toDelete = t.asymmetry & (~((t.Mx == t.My) & (t.period_x == t.period_y) & (t.diam_x == t.diam_y)));
t(toDelete,:) = [];
param_table = t;