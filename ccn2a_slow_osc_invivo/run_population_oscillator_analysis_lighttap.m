%% run_population_oscillator_analysis_lighttap.m
% Population-level oscillator analysis for:
% \\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail
%
% This script:
%   - finds valid fish folders using strict folder logic
%   - loads dF/F and stimulation frames from R
%   - calls oscillator_test_one_recording_v2 for each recording
%   - saves per-recording figures
%   - writes a summary table
%   - generates population-level summary figures
%
% -------------------------------------------------------------------------

clear; clc;

%% ========================== USER SETTINGS ===============================

baseFolder = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail';

cfg.resultsFileName = 'dffs_repact_respcells.mat';
cfg.metaFileName    = 'metadata_multimodal.mat';

cfg.dffVarCandidates = { ...
    'dff_new', ...
    'dff', ...
    'DFF', ...
    'dffMat'};

cfg.fpsVarCandidates = { ...
    'fps', ...
    'fs', ...
    'Fs', ...
    'origFps'};

% stimulus variable names NOW SEARCHED IN R
cfg.stimVarCandidates = { ...
    'stims', ...
    'stim_frames', ...
    'stimOnsetFrames', ...
    'stim_onset_frames', ...
    'stimOnsetsFrames', ...
    'sensoryStimFrames', ...
    'sensorStimFrames', ...
    'lightFrames', ...
    'tapFrames', ...
    'lightTapFrames', ...
    'eventFrames', ...
    'redFrames', ...
    'vibrationFrames', ...
    'combinedFrames'};

% Analysis settings
% Epoch specification
% Options:
%   []                              -> full recording
%   [startFrame endFrame]           -> explicit frame indices
%   struct('unit','min',...)        -> time-based interval, converted per recording using fps
%   struct('unit','sec',...)        -> same idea in seconds
cfg.epochFrames = struct('unit','min', 'start', 2, 'end', 17);
cfg.finiteFracThresh   = 0.95;
cfg.useZScore          = true;
cfg.numPCs             = 10;
cfg.bandHz             = [0.01 0.10];
cfg.welchWinSec        = 120;
cfg.welchOverlapFrac   = 0.5;
cfg.ldsReg             = 1e-3;
cfg.minEigMag          = 0.90;
cfg.maxEigMag          = 1.05;

% Optional ACF oscillation-enrichment filter
cfg.enableAcfPeakFilter = true;
cfg.acfPeakLagRangeSec  = [10 120];
cfg.acfMinProminence    = 0.08;
cfg.minNumAcPeaks       = 2;
cfg.acfMaxLagSec        = 180;

% Mode-membership plot
cfg.makeModeMembershipPlot = true;

% Phase handling
cfg.phaseAmpQuantileMin = 0.20;
cfg.plvThreshold        = 0.70;

% Sliding-window analysis
cfg.winSec             = 600;         % 10 min windows
cfg.stepSec            = 120;         % 2 min steps

% Save / plot
cfg.doPlots            = true;
cfg.saveFigures        = true;
cfg.verbose            = true;

outDir = fullfile(baseFolder, 'oscillator_analysis_results_v2');
figDir = fullfile(outDir, 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% ====================== FIND FISH FOLDERS WITH DATA ====================
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};
resultsFiles = {};
metaFiles    = {};
fishLabel    = {};

for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    rf = fullfile(allFolders{i}, cfg.resultsFileName);
    mf = fullfile(allFolders{i}, cfg.metaFileName);
    if exist(rf,'file')==2 && exist(mf,'file')==2
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
        resultsFiles{end+1} = rf;            %#ok<SAGROW>
        metaFiles{end+1}    = mf;            %#ok<SAGROW>
        [~, fishLabel{end+1}] = fileparts(allFolders{i}); %#ok<SAGROW>
    end
end

fprintf('N fish folders: %d\n', numel(validFolders));

if isempty(validFolders)
    error('No valid folders found.');
end

%% ============================== RUN ====================================

AllRes = struct([]);
summaryRows = [];

for iRec = 1:numel(validFolders)
    fprintf('\n============================================================\n');
    fprintf('Recording %d / %d\n', iRec, numel(validFolders));
    fprintf('Folder: %s\n', validFolders{iRec});

    try
        % ---------- load results file ----------
        R = load(resultsFiles{iRec});
        [dff, fps] = extract_dff_and_fps(R, cfg, resultsFiles{iRec});

        % ---------- stim frames now come from R ----------
        stimFrames = extract_stim_frames(R, cfg);
        stimFrames = unique(round(stimFrames(:)));
        stimFrames = stimFrames(isfinite(stimFrames));
        stimFrames = stimFrames(stimFrames >= 1 & stimFrames <= size(dff,2));
        positions3D = extract_positions_3d(R);   % new helper

        % optional metadata load, retained in case you need metadata later
        M = load(metaFiles{iRec}); %#ok<NASGU>

        % ---------- epoch ----------
        epochFrames = resolve_epoch_spec(cfg.epochFrames, fps, size(dff,2));

        % ---------- run analysis ----------
      OUT = oscillator_test_one_recording_v3( ...
        dff, fps, ...
        'EpochFrames', epochFrames, ...
        'EventFrames', stimFrames, ...
        'FiniteFracThresh', cfg.finiteFracThresh, ...
        'UseZScore', cfg.useZScore, ...
        'NumPCs', cfg.numPCs, ...
        'BandHz', cfg.bandHz, ...
        'WelchWinSec', cfg.welchWinSec, ...
        'WelchOverlapFrac', cfg.welchOverlapFrac, ...
        'LDSReg', cfg.ldsReg, ...
        'MinEigMag', cfg.minEigMag, ...
        'MaxEigMag', cfg.maxEigMag, ...
        'PhaseAmpQuantileMin', cfg.phaseAmpQuantileMin, ...
        'PLVThreshold', cfg.plvThreshold, ...
        'WinSec', cfg.winSec, ...
        'StepSec', cfg.stepSec, ...
        'EnableAcfPeakFilter', cfg.enableAcfPeakFilter, ...
        'AcfPeakLagRangeSec', cfg.acfPeakLagRangeSec, ...
        'AcfMinProminence', cfg.acfMinProminence, ...
        'MinNumAcPeaks', cfg.minNumAcPeaks, ...
        'AcfMaxLagSec', cfg.acfMaxLagSec, ...
        'MakeModeMembershipPlot', cfg.makeModeMembershipPlot, ...
        'ModeSmoothSec', 10, ...
        'ModeDomMarginThr', 0.15,...
        'NeuronPositions', positions3D, ...
        'DoPlots', cfg.doPlots, ...
        'Verbose', cfg.verbose);
            OUT.recordingLabel = fishLabel{iRec};
            OUT.folder = validFolders{iRec};
            OUT.resultsFile = resultsFiles{iRec};
            OUT.metaFile = metaFiles{iRec};

        % ---------- save figures ----------
        if cfg.saveFigures && cfg.doPlots
            save_one_figure_if_present(OUT.fig, 'summary',        figDir, [fishLabel{iRec} '_summary.png']);
            save_one_figure_if_present(OUT.fig, 'phase',          figDir, [fishLabel{iRec} '_phase.png']);
            save_one_figure_if_present(OUT.fig, 'sliding',        figDir, [fishLabel{iRec} '_sliding.png']);
            save_one_figure_if_present(OUT.fig, 'modeMembership', figDir, [fishLabel{iRec} '_modeMembership.png']);
            save_one_figure_if_present(OUT.fig, 'modeDynamics',   figDir, [fishLabel{iRec} '_modeDynamics.png']);
        end

        AllRes = [AllRes; OUT]; %#ok<AGROW>

        % ---------- one summary row per recording ----------
        row = struct();
        row.recordingLabel      = string(fishLabel{iRec});
        row.folder              = string(validFolders{iRec});
        row.nNeurons            = size(dff,1);
        row.nFrames             = size(dff,2);
        row.fps                 = fps;
        row.durationSec         = size(dff,2) / fps;
        row.nStim               = numel(stimFrames);

        row.nKeptNeurons        = OUT.summary.nKeptNeurons;
        row.pcFreqCV            = OUT.summary.pcFreqCV;
        row.meanPLV             = OUT.summary.meanPLV;
        row.nPhaseComponents    = OUT.summary.nPhaseComponents;
        row.nOscillatorPairsLDS = OUT.summary.nOscillatorPairsLDS;
        row.meanSlidingPairs    = mean(OUT.sliding.nOscPairs, 'omitnan');
        row.maxSlidingPairs     = max(OUT.sliding.nOscPairs, [], 'omitnan');
        row.interpretation      = string(OUT.summary.interpretation);

        row.pc1_domFreqHz       = safe_idx(OUT.summary.pcDomFreqHz, 1);
        row.pc2_domFreqHz       = safe_idx(OUT.summary.pcDomFreqHz, 2);
        row.pc3_domFreqHz       = safe_idx(OUT.summary.pcDomFreqHz, 3);
        row.pc4_domFreqHz       = safe_idx(OUT.summary.pcDomFreqHz, 4);
        row.pc5_domFreqHz       = safe_idx(OUT.summary.pcDomFreqHz, 5);
        row.pc6_domFreqHz       = safe_idx(OUT.summary.pcDomFreqHz, 6);

        row.expVar_PC1          = safe_idx(OUT.pca.explained, 1);
        row.expVar_PC2          = safe_idx(OUT.pca.explained, 2);
        row.expVar_PC3          = safe_idx(OUT.pca.explained, 3);
        row.expVar_PC4          = safe_idx(OUT.pca.explained, 4);
        row.expVar_PC5          = safe_idx(OUT.pca.explained, 5);
        row.expVar_PC6          = safe_idx(OUT.pca.explained, 6);

        row.nAcfFilteredOut     = sum(~OUT.filter.keepFinal);
        row.nAcfKept            = sum(OUT.filter.keepFinal);

        row.mode1_n             = sum(OUT.mode.modeID == 1);
        row.mode2_n             = sum(OUT.mode.modeID == 2);
        row.mode3_n             = sum(OUT.mode.modeID == 3);
        row.mode4_n             = sum(OUT.mode.modeID == 4);
        row.mode5_n             = sum(OUT.mode.modeID == 5);
        row.mode6_n             = sum(OUT.mode.modeID == 6);

        row.nOccupiedModes      = nnz(OUT.summary.modeCount > 0);

        summaryRows = [summaryRows; row]; %#ok<AGROW>

    catch ME
        warning('Failed on folder:\n%s\n\n%s', validFolders{iRec}, ...
            getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

%% ============================== SAVE ===================================
if isempty(AllRes)
    error('No recordings completed successfully.');
end

summaryTable = struct2table(summaryRows);

writetable(summaryTable, fullfile(outDir, 'summary_table.csv'));
save(fullfile(outDir, 'summary_struct.mat'), 'AllRes', 'summaryTable', 'cfg', ...
    'validFolders', 'resultsFiles', 'metaFiles', 'fishLabel', '-v7.3');

disp(summaryTable);

%% ===================== POPULATION SUMMARY PLOTS ========================
popFigs = make_population_summary_plots(summaryTable, AllRes);

if cfg.saveFigures
    save_population_fig_if_valid(popFigs, 'overview',    figDir, 'population_overview.png');
    save_population_fig_if_valid(popFigs, 'frequencies', figDir, 'population_frequencies.png');
    save_population_fig_if_valid(popFigs, 'variance',    figDir, 'population_explained_variance.png');
    save_population_fig_if_valid(popFigs, 'modes',       figDir, 'population_modes.png');
end

%% =========================== LOCAL FUNCTIONS ===========================

function save_one_figure_if_present(figStruct, fieldName, figDir, fileName)
if isfield(figStruct, fieldName)
    figObj = figStruct.(fieldName);
    if ~isempty(figObj) && isgraphics(figObj)
        exportgraphics(figObj, fullfile(figDir, fileName), 'Resolution', 200);
    end
end
end

function save_population_fig_if_valid(figStruct, fieldName, figDir, fileName)
if isfield(figStruct, fieldName)
    figObj = figStruct.(fieldName);
    if ~isempty(figObj) && isgraphics(figObj)
        exportgraphics(figObj, fullfile(figDir, fileName), 'Resolution', 200);
    end
end
end

function val = safe_idx(x, idx)
if numel(x) >= idx
    val = x(idx);
else
    val = NaN;
end
end

function [dff, fps] = extract_dff_and_fps(R, cfg, filePath)
dff = [];
fps = [];

vars = fieldnames(R);

for i = 1:numel(cfg.dffVarCandidates)
    vn = cfg.dffVarCandidates{i};
    if isfield(R, vn)
        cand = R.(vn);
        if isnumeric(cand) && ismatrix(cand) && size(cand,1) > 1 && size(cand,2) > 10
            dff = cand;
            break;
        end
    end
end

if isempty(dff)
    for i = 1:numel(vars)
        x = R.(vars{i});
        if isstruct(x)
            nested = {'dff_new','dff','DFF','DV_DFFmovwindow'};
            for j = 1:numel(nested)
                if isfield(x, nested{j})
                    cand = x.(nested{j});
                    if isnumeric(cand) && ismatrix(cand) && size(cand,1) > 1 && size(cand,2) > 10
                        dff = cand;
                        break;
                    end
                end
            end
        end
        if ~isempty(dff), break; end
    end
end

if isempty(dff)
    error('Could not find dF/F matrix in %s', filePath);
end

for i = 1:numel(cfg.fpsVarCandidates)
    vn = cfg.fpsVarCandidates{i};
    if isfield(R, vn) && isnumeric(R.(vn)) && isscalar(R.(vn))
        fps = double(R.(vn));
        break;
    end
end

if isempty(fps)
    for i = 1:numel(vars)
        x = R.(vars{i});
        if isstruct(x)
            for j = 1:numel(cfg.fpsVarCandidates)
                vn = cfg.fpsVarCandidates{j};
                if isfield(x, vn) && isnumeric(x.(vn)) && isscalar(x.(vn))
                    fps = double(x.(vn));
                    break;
                end
            end
        end
        if ~isempty(fps), break; end
    end
end

if isempty(fps)
    error('Could not find fps in %s', filePath);
end
end

function stimFrames = extract_stim_frames(S, cfg)
stimFrames = [];
vars = fieldnames(S);

for i = 1:numel(cfg.stimVarCandidates)
    vn = cfg.stimVarCandidates{i};
    if isfield(S, vn)
        cand = S.(vn);
        if isnumeric(cand) && isvector(cand)
            stimFrames = cand(:);
            return;
        end
    end
end

for i = 1:numel(vars)
    x = S.(vars{i});
    if isstruct(x)
        for j = 1:numel(cfg.stimVarCandidates)
            vn = cfg.stimVarCandidates{j};
            if isfield(x, vn)
                cand = x.(vn);
                if isnumeric(cand) && isvector(cand)
                    stimFrames = cand(:);
                    return;
                end
            end
        end
    end
end
end

function figs = make_population_summary_plots(summaryTable, AllRes)
figs = struct();

% ===================== OVERVIEW =====================
figs.overview = figure('Color','w','Position',[100 100 1500 900], 'Name','Population overview');
tlo = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile(tlo);
histogram(ax1, summaryTable.nOscillatorPairsLDS, 'BinMethod', 'integers');
xlabel(ax1, '# LDS oscillator pairs');
ylabel(ax1, '# recordings');
title(ax1, 'Whole-recording oscillator pair counts');

ax2 = nexttile(tlo);
scatter(ax2, summaryTable.pcFreqCV, summaryTable.meanPLV, 50, 'filled');
xlabel(ax2, 'PC frequency CV');
ylabel(ax2, 'mean PLV');
title(ax2, 'Timescale spread vs phase locking');
grid(ax2, 'on');

ax3 = nexttile(tlo);
boxchart(ax3, summaryTable.nOccupiedModes);
ylabel(ax3, '# occupied modes');
title(ax3, 'Occupied modes across recordings');

ax4 = nexttile(tlo);
histogram(ax4, summaryTable.maxSlidingPairs, 'BinMethod', 'integers');
xlabel(ax4, 'max sliding-window oscillator pairs');
ylabel(ax4, '# recordings');
title(ax4, 'Transient oscillator evidence');

% ===================== FREQUENCIES =====================
figs.frequencies = figure('Color','w','Position',[120 120 1500 900], 'Name','Population frequencies');
tlo2 = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile(tlo2);
boxchart(ax1, [summaryTable.pc1_domFreqHz, summaryTable.pc2_domFreqHz, summaryTable.pc3_domFreqHz, ...
               summaryTable.pc4_domFreqHz, summaryTable.pc5_domFreqHz, summaryTable.pc6_domFreqHz]);
xlabel(ax1, 'PC');
ylabel(ax1, 'dominant frequency (Hz)');
title(ax1, 'PC dominant frequencies across recordings');

ax2 = nexttile(tlo2);
histogram(ax2, summaryTable.pcFreqCV);
xlabel(ax2, 'PC frequency CV');
ylabel(ax2, '# recordings');
title(ax2, 'Spread of latent timescales');

ax3 = nexttile(tlo2, [1 2]);
hold(ax3, 'on');
for i = 1:height(summaryTable)
    vals = [summaryTable.pc1_domFreqHz(i), summaryTable.pc2_domFreqHz(i), summaryTable.pc3_domFreqHz(i), ...
            summaryTable.pc4_domFreqHz(i), summaryTable.pc5_domFreqHz(i), summaryTable.pc6_domFreqHz(i)];
    plot(ax3, 1:6, vals, '-o', 'LineWidth', 1);
end
xlabel(ax3, 'PC');
ylabel(ax3, 'dominant frequency (Hz)');
title(ax3, 'PC dominant frequencies by recording');
grid(ax3, 'on');

% ===================== EXPLAINED VARIANCE =====================
figs.variance = figure('Color','w','Position',[140 140 1500 900], 'Name','Population explained variance');
tlo3 = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

expMat = [summaryTable.expVar_PC1, summaryTable.expVar_PC2, summaryTable.expVar_PC3, ...
          summaryTable.expVar_PC4, summaryTable.expVar_PC5, summaryTable.expVar_PC6];

ax1 = nexttile(tlo3);
boxchart(ax1, expMat);
xlabel(ax1, 'PC');
ylabel(ax1, '% variance explained');
title(ax1, 'Explained variance across recordings');

ax2 = nexttile(tlo3);
bar(ax2, mean(expMat,1,'omitnan'));
xlabel(ax2, 'PC');
ylabel(ax2, 'mean % variance explained');
title(ax2, 'Mean explained variance');

ax3 = nexttile(tlo3, [1 2]);
imagesc(ax3, expMat);
xlabel(ax3, 'PC');
ylabel(ax3, 'recording');
title(ax3, 'Explained variance heatmap');
colorbar(ax3);

% ===================== MODES =====================
figs.modes = figure('Color','w','Position',[160 160 1500 900], 'Name','Population mode summaries');
tlo4 = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

modeMat = [summaryTable.mode1_n, summaryTable.mode2_n, summaryTable.mode3_n, ...
           summaryTable.mode4_n, summaryTable.mode5_n, summaryTable.mode6_n];

ax1 = nexttile(tlo4);
boxchart(ax1, modeMat);
xlabel(ax1, 'mode');
ylabel(ax1, '# neurons');
title(ax1, 'Mode sizes across recordings');

ax2 = nexttile(tlo4);
imagesc(ax2, modeMat);
xlabel(ax2, 'mode');
ylabel(ax2, 'recording');
title(ax2, 'Mode-size heatmap');
colorbar(ax2);

ax3 = nexttile(tlo4);
histogram(ax3, summaryTable.nOccupiedModes, 'BinMethod', 'integers');
xlabel(ax3, '# occupied modes');
ylabel(ax3, '# recordings');
title(ax3, 'Distribution of occupied modes');

ax4 = nexttile(tlo4);
scatter(ax4, summaryTable.nAcfKept, summaryTable.nOccupiedModes, 50, 'filled');
xlabel(ax4, '# neurons kept after ACF filter');
ylabel(ax4, '# occupied modes');
title(ax4, 'ACF-kept neurons vs occupied modes');
grid(ax4, 'on');

end

function epochFrames = resolve_epoch_spec(epochSpec, fps, nFramesTotal)
% RESOLVE_EPOCH_SPEC
% Convert cfg.epochFrames into [startFrame endFrame].
%
% Accepted inputs:
%   []                              -> full recording
%   [startFrame endFrame]           -> explicit frame indices
%   struct('unit','min',...)        -> fields: start, end
%   struct('unit','sec',...)        -> fields: start, end

    if isempty(epochSpec)
        epochFrames = [1 nFramesTotal];
        return;
    end

    % Case 1: explicit frame indices
    if isnumeric(epochSpec) && numel(epochSpec) == 2
        epochFrames = round(epochSpec(:))';
        epochFrames(1) = max(1, epochFrames(1));
        epochFrames(2) = min(nFramesTotal, epochFrames(2));

        if epochFrames(2) <= epochFrames(1)
            error('Resolved epoch frame interval is invalid.');
        end
        return;
    end

    % Case 2: time-based struct
    if isstruct(epochSpec)
        if ~isfield(epochSpec, 'unit') || ~isfield(epochSpec, 'start') || ~isfield(epochSpec, 'end')
            error('epochSpec struct must contain fields: unit, start, end.');
        end

        unitStr = lower(string(epochSpec.unit));
        tStart = double(epochSpec.start);
        tEnd   = double(epochSpec.end);

        switch unitStr
            case "min"
                startFrame = floor(tStart * 60 * fps) + 1;
                endFrame   = floor(tEnd   * 60 * fps);

            case "sec"
                startFrame = floor(tStart * fps) + 1;
                endFrame   = floor(tEnd   * fps);

            otherwise
                error('epochSpec.unit must be ''min'' or ''sec''.');
        end

        epochFrames = [max(1, startFrame), min(nFramesTotal, endFrame)];

        if epochFrames(2) <= epochFrames(1)
            error('Resolved epoch frame interval is invalid.');
        end
        return;
    end

    error('cfg.epochFrames must be [], a 2-element numeric vector, or a struct with unit/start/end.');
end

function positions3D = extract_positions_3d(S)
positions3D = [];

% direct fields
candNames = {'position','positions','pos','pos3D','xyz'};
for i = 1:numel(candNames)
    if isfield(S, candNames{i})
        x = S.(candNames{i});
        if isnumeric(x) && ismatrix(x) && (size(x,2)==3 || size(x,2)==5)
            positions3D = x(:,1:3);
            break;
        end
    end
end

% nested struct fields
if isempty(positions3D)
    vars = fieldnames(S);
    for i = 1:numel(vars)
        x = S.(vars{i});
        if isstruct(x)
            nestedNames = {'position','positions3D','positions','xyz'};
            for j = 1:numel(nestedNames)
                if isfield(x, nestedNames{j})
                    p = x.(nestedNames{j});
                    if isnumeric(p) && ismatrix(p) && (size(p,2)==3 || size(p,2)==5)
                        positions3D = p(:,1:3);
                        break;
                    end
                end
            end
        end
        if ~isempty(positions3D)
            break;
        end
    end
end

if isempty(positions3D)
    positions3D = [];
    return;
end

% if only XY, append Z=0
if size(positions3D,2) == 2
    positions3D = [positions3D, zeros(size(positions3D,1),1)];
end
end