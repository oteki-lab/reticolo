period_x = [0; 50];
period_y = [0; 15];
score = [44; 68];
hp_steps = table(period_x,period_y,score);
%steps = [
%   [50, 15, 68], %初期点
%   [0, 0, 44],   %初期点
%   [25, 15, 58], %実験1回目のデータ
%   [50, 0, 85],  %2回目
%   [45, 0, 94],  %3回目...
%   [25, 6, 50],
%   [50, 6, 90],
%   [0, 15, 46],
%];
save('hp_steps.mat','hp_steps');
writetable(hp_steps, "hp_steps.csv", 'Delimiter',',');