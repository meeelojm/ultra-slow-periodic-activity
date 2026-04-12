%% --------- Define your files (fish × {high,mid,low}) ----------
files = struct([]);

% F012
files(1).name = "F012";
files(1).high = "Z:\temp\Emilian temp\ccn2a\high2low\f012\high\ethogram_20250828_065717.mat";
files(1).mid  = "Z:\temp\Emilian temp\ccn2a\high2low\f012\mid\ethogram_20250828_132443.mat";
files(1).low  = "Z:\temp\Emilian temp\ccn2a\high2low\f012\low\ethogram_20250828_204355.mat";

% F010
files(2).name = "F010";
files(2).high = "Z:\temp\Emilian temp\ccn2a\high2low\f010\highact\suite2p\ethogram_20250820_110948_2.mat";
files(2).mid  = "Z:\temp\Emilian temp\ccn2a\high2low\f010\midact\suite2p\ethogram_20250821_132810.mat";
files(2).low  = "Z:\temp\Emilian temp\ccn2a\high2low\f010\lowact\suite2p\ethogram_20250821_115816.mat";

% F011
files(3).name = "F011";
files(3).high = "Z:\temp\Emilian temp\ccn2a\high2low\f011\highact\suite2p\ethogram_20250820_195637.mat";
files(3).mid  = "Z:\temp\Emilian temp\ccn2a\high2low\f011\midact\suite2p\ethogram_20250820_103528.mat";
files(3).low  = "Z:\temp\Emilian temp\ccn2a\high2low\f011\lowact\suite2p\ethogram_20250820_104458.mat";

nA = numel(files);

%% --------- Load each (time, angle) ----------
H = cell(nA,1); M = cell(nA,1); L = cell(nA,1);
tH = cell(nA,1); tM = cell(nA,1); tL = cell(nA,1);

for a = 1:nA
    [tH{a}, H{a}] = load_eth(files(a).high);
    [tM{a}, M{a}] = load_eth(files(a).mid);
    [tL{a}, L{a}] = load_eth(files(a).low);
end

%% --------- Build a single common time base (intersection) ----------
starts = cellfun(@(x) x(1), [tH; tM; tL]);
ends   = cellfun(@(x) x(end), [tH; tM; tL]);
dt_all = cellfun(@(x) median(diff(x)), [tH; tM; tL]);

t0 = max(starts);
t1 = min(ends);
assert(t1 > t0, 'No time overlap across all traces.');
dt_common = median(dt_all);                     % could use min(dt_all) if you prefer finer grid
t_common  = (t0:dt_common:t1).';

% Resample & stack columns = fish


for a = 1:nA
    high_mat(:,a) = resamp_to(tH{a}, H{a}, t_common);
    mid_mat(:,a)  = resamp_to(tM{a}, M{a}, t_common);
    low_mat(:,a)  = resamp_to(tL{a}, L{a}, t_common);
end

%% --------- Quantify across conditions ----------
params = struct('fs', 1/median(diff(t_common)), ...
                'smooth_sec', 0.10, ...
                'thr_k', 4, ...
                'min_bout_sec', 0.10);
OUT = quantify_tail_by_condition(t_common, high_mat, mid_mat, low_mat, params);

