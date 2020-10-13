period_x = [50; 0];
period_y = [15; 0];
score = [68; 44];
steps = table(period_x,period_y,score);
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
save('step_list.mat','steps');