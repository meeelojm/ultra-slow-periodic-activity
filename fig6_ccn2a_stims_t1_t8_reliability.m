%% --- 
% Define the base folder path
baseFolder = 'Z:\mh-kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';

% Get a list of all subfolders (including the base folder) using genpath
allFolders = strsplit(genpath(baseFolder), pathsep);

% Initialize a cell array to store folders that contain both files
validFolders = {};

% Loop over each folder and check for the required files
for i = 1:length(allFolders)
    currentFolder = allFolders{i};
    if isempty(currentFolder)
        continue; % Skip empty entries
    end
    
    % Construct the full file paths for the two files
    resultsFile = fullfile(currentFolder, 'dffs_repact_respcells.mat');
    metadataFile = fullfile(currentFolder, 'metadata_multimodal.mat');
    
    % Check if both files exist (exist returns 2 for files)
    if exist(resultsFile, 'file') == 2 %&& exist(metadataFile, 'file') == 2
        validFolders{end+1} = currentFolder;  %#ok<SAGROW>
    end
end

% Display the folders that contain both files
disp('Folders containing both results.mat and metadata_multimodal.mat:');
for i = 1:length(validFolders)
    disp(validFolders{i});
end

%%

% Preallocate cell arrays for each condition (for "p" and "n" groups)
all_plight = cell(length(validFolders), 1);
all_ptap   = cell(length(validFolders), 1);
all_plt    = cell(length(validFolders), 1);

all_nlight = cell(length(validFolders), 1);
all_ntap   = cell(length(validFolders), 1);
all_nlt    = cell(length(validFolders), 1);

all_nsize  = cell(length(validFolders), 1);

for i = 1:length(validFolders)
    dataFolder = validFolders{i}; 
    % Load only the required variables from the file.
    temp = load(fullfile(dataFolder, "dffs_repact_respcells.mat"), ...
        'sorting_mean_plight_dstd', 'sorting_mean_ptap_dstd', 'sorting_mean_plt_dstd', ...
        'sorting_mean_nlight_dstd', 'sorting_mean_ntap_dstd', 'sorting_mean_nlt_dstd', 'dff');
    
    % For the "p" group, transpose to get time x cells, then take the first 600 rows.
    all_plight{i} = temp.sorting_mean_plight_dstd(1:500, :)';
    all_ptap{i}   = temp.sorting_mean_ptap_dstd(1:500, :)';
    all_plt{i}    = temp.sorting_mean_plt_dstd(1:500, :)';
    
    % For the "n" group.
    all_nlight{i} = temp.sorting_mean_nlight_dstd(1:500, :)';
    all_ntap{i}   = temp.sorting_mean_ntap_dstd(1:500, :)';
    all_nlt{i}    = temp.sorting_mean_nlt_dstd(1:500, :)';
    
    % Store the number of cells from dff.
    all_nsize{i} = size(temp.dff, 1);
end

% Concatenate data from all folders (this is much faster than growing arrays in the loop)
sm_plight = vertcat(all_plight{:});
sm_ptap   = vertcat(all_ptap{:});
sm_plt    = vertcat(all_plt{:});

sm_nlight = vertcat(all_nlight{:});
sm_ntap   = vertcat(all_ntap{:});
sm_nlt    = vertcat(all_nlt{:});

% Convert the cell array of sizes to a numeric vector.
nsize = cell2mat(all_nsize);


%%
% Preallocate cell arrays for each condition (S1 and S8)
all_light_s1 = cell(length(validFolders),1);
all_tap_s1   = cell(length(validFolders),1);
all_lt_s1    = cell(length(validFolders),1);

all_light_s8 = cell(length(validFolders),1);
all_tap_s8   = cell(length(validFolders),1);
all_lt_s8    = cell(length(validFolders),1);

for i = 1:length(validFolders)
    dataFolder = validFolders{i}; 
    % Load only the necessary variables from the file
    temp = load(fullfile(dataFolder, "dffs_repact_respcells.mat"), ...
        'sorting_mean_s1_alllight_dstd','sorting_mean_s1_alltap_dstd','sorting_mean_s1_alllt_dstd',...
        'sorting_mean_s8_alllight_dstd','sorting_mean_s8_alltap_dstd','sorting_mean_s8_alllt_dstd');
    
    % For S1: Transpose to get time x cells, then take the first 600 rows
    all_light_s1{i} = temp.sorting_mean_s1_alllight_dstd(1:500,:)';
    all_tap_s1{i}   = temp.sorting_mean_s1_alltap_dstd(1:500,:)';
    all_lt_s1{i}    = temp.sorting_mean_s1_alllt_dstd(1:500,:)';
    
    % For S8:
    all_light_s8{i} = temp.sorting_mean_s8_alllight_dstd(1:500,:)';
    all_tap_s8{i}   = temp.sorting_mean_s8_alltap_dstd(1:500,:)';
    all_lt_s8{i}    = temp.sorting_mean_s8_alllt_dstd(1:500,:)';
end
assert(all(cellfun(@(a,b) size(a,1)==size(b,1), all_light_s1, all_light_s8)), ...
    'Mismatch in neuron counts between S1 and S8 for some fish.');
% Use the S1 matrices to define nsize (rows = neurons)
nsize = cellfun(@(M) size(M,1), all_light_s1);   % 1 x nFish (or column if you prefer)
nsize = nsize(:);  
%%
% Concatenate data from all folders (this is much faster than appending in the loop)
sm_alllight_s1 = vertcat(all_light_s1{:});
sm_alltap_s1   = vertcat(all_tap_s1{:});
sm_alllt_s1    = vertcat(all_lt_s1{:});

sm_alllight_s8 = vertcat(all_light_s8{:});
sm_alltap_s8   = vertcat(all_tap_s8{:});
sm_alllt_s8    = vertcat(all_lt_s8{:});

% ---- Light
[rel_light, rel_light_null, perFish_light] = split_reliability(sm_alllight_s1, sm_alllight_s8, nsize);

% ---- Tap
[rel_tap, rel_tap_null, perFish_tap] = split_reliability(sm_alltap_s1, sm_alltap_s8, nsize);

% ---- Light+Tap
[rel_lt, rel_lt_null, perFish_lt] = split_reliability(sm_alllt_s1, sm_alllt_s8, nsize);

% ---- Quick summaries
fprintf('\n=== S1↔S8 reliability (median [IQR]) ===\n');
summ = @(x) sprintf('%.3f  [%.3f–%.3f]', median(x,'omitnan'), quantile(x,0.25), quantile(x,0.75));
fprintf('Light   : %s  | null: %s\n', summ(rel_light), summ(rel_light_null));
fprintf('Tap     : %s  | null: %s\n', summ(rel_tap),   summ(rel_tap_null));
fprintf('Light+Tap: %s | null: %s\n', summ(rel_lt),    summ(rel_lt_null));

% ---- Per-fish medians (handy for stats/plots)
perFishMeds_light = cellfun(@(v) median(v,'omitnan'), perFish_light);
perFishMeds_tap   = cellfun(@(v) median(v,'omitnan'), perFish_tap);
perFishMeds_lt    = cellfun(@(v) median(v,'omitnan'), perFish_lt);

% ---- Plots
figure('Color','w','Name','S1↔S8 reliability'); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; histogram(rel_light, -1:0.05:1); hold on; histogram(rel_light_null, -1:0.05:1);
title('Light'); xlabel('r'); ylabel('# neurons'); legend('real','shuffle'); box on;

nexttile; histogram(rel_tap, -1:0.05:1); hold on; histogram(rel_tap_null, -1:0.05:1);
title('Tap'); xlabel('r'); ylabel('# neurons'); legend('real','shuffle'); box on;

nexttile; histogram(rel_lt, -1:0.05:1); hold on; histogram(rel_lt_null, -1:0.05:1);
title('Light+Tap'); xlabel('r'); ylabel('# neurons'); legend('real','shuffle'); box on;

% Per-fish box/strip (medians) — absolute
vL = abs(perFishMeds_light);
vT = abs(perFishMeds_tap);
vLT = abs(perFishMeds_lt);

nexttile; boxchart(ones(numel(vL),1), vL,'BoxFaceAlpha',0.2); hold on;
scatter(ones(numel(vL),1)+0.05*randn(numel(vL),1), vL, 18, 'filled','MarkerFaceAlpha',0.7);
xticks([]); ylabel('median |r|'); title('Per-fish Light');

nexttile; boxchart(ones(numel(vT),1), vT,'BoxFaceAlpha',0.2); hold on;
scatter(ones(numel(vT),1)+0.05*randn(numel(vT),1), vT, 18, 'filled','MarkerFaceAlpha',0.7);
xticks([]); ylabel('median |r|'); title('Per-fish Tap');

nexttile; boxchart(ones(numel(vLT),1), vLT,'BoxFaceAlpha',0.2); hold on;
scatter(ones(numel(vLT),1)+0.05*randn(numel(vLT),1), vLT, 18, 'filled','MarkerFaceAlpha',0.7);
xticks([]); ylabel('median |r|'); title('Per-fish Light+Tap');


%%  Mapping (Sorting) Using mapTmap
ops.nC = 5;
ops.iPC = 1:10;
ops.isort = [];
ops.useGPU = 1;
ops.upsamp = 100;
ops.sigUp = 1;

% Process S1 data for each condition
[isort1_alllight_s1, ~, ~, ~, iclust] = mapTmap(sm_alllight_s1, ops);
sm_alllight_s1s1 = sm_alllight_s1(isort1_alllight_s1, :);

[isort1_alltap_s1, ~, ~, ~, iclust] = mapTmap(sm_alltap_s1, ops);
sm_alltap_s1s1 = sm_alltap_s1(isort1_alltap_s1, :);

[isort1_alllt_s1, ~, ~, ~, iclust] = mapTmap(sm_alllt_s1, ops);
sm_alllt_s1s1 = sm_alllt_s1(isort1_alllt_s1, :);

% Process S8 data for each condition
[isort1_alllight_s8, ~, ~, ~, iclust] = mapTmap(sm_alllight_s8, ops);
sm_alllight_s8s8 = sm_alllight_s8(isort1_alllight_s8, :);
sm_alllight_s8s1 = sm_alllight_s8(isort1_alllight_s1, :);

[isort1_alltap_s8, ~, ~, ~, iclust] = mapTmap(sm_alltap_s8, ops);
sm_alltap_s8s8 = sm_alltap_s8(isort1_alltap_s8, :);
sm_alltap_s8s1 = sm_alltap_s8(isort1_alltap_s1, :);

[isort1_alllt_s8, ~, ~, ~, iclust] = mapTmap(sm_alllt_s8, ops);
sm_alllt_s8s8 = sm_alllt_s8(isort1_alllt_s8, :);
sm_alllt_s8s1 = sm_alllt_s8(isort1_alllt_s1, :);

% Define condition names for first trial (S1) and last trial (S8)
condNamesS1S1 = {'modulated light', 'modulated tap', 'modulated light + tap'};
dataS1S1 = {sm_alllight_s1s1, sm_alltap_s1s1, sm_alllt_s1s1};

condNamesS8S8 = {'modulated light', 'modulated tap', 'modulated light + tap'};
dataS8S8 = {sm_alllight_s8s8, sm_alltap_s8s8, sm_alllt_s8s8};

condNamesS8S1 = {'modulated light', 'modulated tap', 'modulated light + tap'};
dataS8S1 = {sm_alllight_s8s1, sm_alltap_s8s1, sm_alllt_s8s1};

% Figure for First Trial (S1)
fig1 = figure('Renderer','painters','Position',[820,420,820,420]);
t1 = tiledlayout(fig1, 1, 10, 'TileSpacing','compact','Padding','compact');
createSubplots(fig1, t1, dataS1S1, condNamesS1S1, 'first trial');

% Figure for Last Trial (S8)
fig2 = figure('Renderer','painters','Position',[820,420,820,420]);
t2 = tiledlayout(fig2, 1, 10, 'TileSpacing','compact','Padding','compact');
createSubplots(fig2, t2, dataS8S8, condNamesS8S8, 'last trial');

% Figure for Last Trial (S8)
fig3 = figure('Renderer','painters','Position',[820,420,820,420]);
t3 = tiledlayout(fig3, 1, 10, 'TileSpacing','compact','Padding','compact');
createSubplots(fig3, t3, dataS8S1, condNamesS8S1, 'last trial (sorted by first)');

%% This is the one wwithh all accumulated trials
% Preallocate cell arrays for each condition (S1 and S8)
all_light = cell(length(validFolders),1);
all_tap   = cell(length(validFolders),1);
all_lt    = cell(length(validFolders),1);


for i = 1:length(validFolders)
    dataFolder = validFolders{i}; 
    % Load only the necessary variables from the file
    temp = load(fullfile(dataFolder, "dffs_repact_respcells.mat"), ...
        'sorting_mean_alllight_dstd','sorting_mean_alltap_dstd','sorting_mean_alllt_dstd');
        
    % For all: Transpose to get time x cells, then take the first 600 rows
    all_light{i} = temp.sorting_mean_alllight_dstd(1:500,:)';
    all_tap{i}   = temp.sorting_mean_alltap_dstd(1:500,:)';
    all_lt{i}    = temp.sorting_mean_alllt_dstd(1:500,:)';
    
end

% Concatenate data from all folders (this is much faster than appending in the loop)
sm_alllight = vertcat(all_light{:});
sm_alltap   = vertcat(all_tap{:});
sm_alllt    = vertcat(all_lt{:});

%  Mapping (Sorting) Using mapTmap
ops.nC = 5;
ops.iPC = 1:10;
ops.isort = [];
ops.useGPU = 1;
ops.upsamp = 100;
ops.sigUp = 1;

% Process S1 data for each condition
[isort1_alllight, ~, ~, ~, iclust] = mapTmap(sm_alllight, ops);
sm_alllight = sm_alllight(isort1_alllight, :);

[isort1_alltap, ~, ~, ~, iclust] = mapTmap(sm_alltap, ops);
sm_alltap = sm_alltap(isort1_alltap, :);

[isort1_alllt, ~, ~, ~, iclust] = mapTmap(sm_alllt, ops);
sm_alllt = sm_alllt_s1(isort1_alllt, :);

% Define condition names for first trial (S1) and last trial (S8)
condNames = {'modulated light', 'modulated tap', 'modulated light + tap'};
data = {sm_alllight, sm_alltap, sm_alllt};

% Figure for First Trial (S1)
fig1 = figure('Renderer','painters','Position',[820,420,820,420]);
t1 = tiledlayout(fig1, 1, 10, 'TileSpacing','compact','Padding','compact');
createSubplots(fig1, t1, data, condNames, 'all trials');

%%

sm_alllight_s1all = sm_alllight_s1(isort1_alllight, :);
sm_alltap_s1all = sm_alltap_s1(isort1_alltap, :);
sm_alllt_s1all = sm_alllt_s1(isort1_alllt, :);

sm_alllight_s8all = sm_alllight_s8(isort1_alllight, :);
sm_alltap_s8all = sm_alltap_s8(isort1_alltap, :);
sm_alllt_s8all = sm_alllt_s8(isort1_alllt, :);

% Define condition names for first trial (S1) and last trial (S8)
condNames = {'modulated light', 'modulated tap', 'modulated light + tap'};
data = {sm_alllight_s1all, sm_alltap_s1all, sm_alllt_s1all};

% Figure for First Trial (S1)
fig1 = figure('Renderer','painters','Position',[820,420,820,420]);
t1 = tiledlayout(fig1, 1, 10, 'TileSpacing','compact','Padding','compact');
createSubplots(fig1, t1, data, condNames, 'first trials (sorted by avg)');


% Define condition names for first trial (S1) and last trial (S8)
condNames = {'modulated light', 'modulated tap', 'modulated light + tap'};
data = {sm_alllight_s8all, sm_alltap_s8all, sm_alllt_s8all};

% Figure for First Trial (S1)
fig1 = figure('Renderer','painters','Position',[820,420,820,420]);
t1 = tiledlayout(fig1, 1, 10, 'TileSpacing','compact','Padding','compact');
createSubplots(fig1, t1, data, condNames, 'last trials (sorted by avg)');

%% cross correlation values

for i = 1:length(validFolders)
    dataFolder = validFolders{i}; 
    % Load only the necessary variables from the file
    temp = load(fullfile(dataFolder, "outputs_rastermapmat"), ...
        'sorting_mean_alllight_dstd','sorting_mean_alltap_dstd','sorting_mean_alllt_dstd');
        
    % For all: Transpose to get time x cells, then take the first 600 rows
    all_light{i} = temp.sorting_mean_alllight_dstd(1:600,:)';
    all_tap{i}   = temp.sorting_mean_alltap_dstd(1:600,:)';
    all_lt{i}    = temp.sorting_mean_alllt_dstd(1:600,:)';
    
end

figure
plot(mean(outputs.crossres_ongoing.xcfs,1))
hold on
plot(mean(outputs.crossres_light.xcfs,1))
hold on
plot(mean(outputs.crossres_tap.xcfs,1))
hold on
plot(mean(outputs.crossres_lt.xcfs,1))
hold on
plot(mean(outputs.crossres_outgoing.xcfs,1))
hold on
plot(mean(outputs.crossres_endgoing.xcfs,1))
legend('ongoing','light','tap','l+t','outgoing','endgoing','Location','best');
ylabel('cross correlation');
xline(11,'k')



