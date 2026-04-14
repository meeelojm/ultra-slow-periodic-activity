function OUT = compare_numPCs_group_summary(baseFolder, varargin)
%COMPARE_NUMPCS_GROUP_SUMMARY
% Group-level sweep over latent dimensionality for oscillator_test_one_recording_v2.
%
% Required:
%   baseFolder : root folder containing recording subfolders
%
% This function:
%   1) finds valid folders containing both:
%        - dffs_repact_respcells.mat
%        - metadata_multimodal.mat
%   2) runs oscillator_test_one_recording_v2 for numPCs in [3 4 5 6 8] by default
%   3) computes compact group summaries to help choose dimensionality
%
% Main outputs per numPCs:
%   - explained variance retained
%   - occupied modes
%   - mean PLV
%   - nPhaseComponents
%   - oscillator pairs (whole-recording and sliding-window)
%   - geometric separation of mode-dominated points in latent space
%   - stability of neuron assignments vs previous numPCs
%
% Example:
%   OUT = compare_numPCs_group_summary( ...
%       '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\light tap free tail', ...
%       'NumPCsList', [3 4 5 6 8], ...
%       'EpochSpec', struct('unit','min','start',2,'end',27), ...
%       'EnableAcfPeakFilter', true, ...
%       'AcfPeakLagRangeSec', [10 120], ...
%       'AcfMinProminence', 0.08, ...
%       'MinNumAcPeaks', 2, ...
%       'Verbose', true);
%
% -------------------------------------------------------------------------

%% ============================= PARSE ===================================
ip = inputParser;
addRequired(ip, 'baseFolder', @(x) ischar(x) || isstring(x));

addParameter(ip, 'ResultsFileName', 'dffs_repact_respcells.mat', @(x) ischar(x) || isstring(x));
addParameter(ip, 'MetaFileName', 'metadata_multimodal.mat', @(x) ischar(x) || isstring(x));

addParameter(ip, 'DffVarCandidates', {'dff_new','dff','DFF','dffMat'}, @iscell);
addParameter(ip, 'FpsVarCandidates', {'fps','fs','Fs','origFps'}, @iscell);
addParameter(ip, 'StimVarCandidates', { ...
    'stims','stim_frames','stimOnsetFrames','stim_onset_frames', ...
    'stimOnsetsFrames','sensoryStimFrames','sensorStimFrames', ...
    'lightFrames','tapFrames','lightTapFrames','eventFrames', ...
    'redFrames','vibrationFrames','combinedFrames'}, @iscell);

addParameter(ip, 'NumPCsList', [3 4 5 6 8], @(x) isnumeric(x) && isvector(x) && all(x>=2));
addParameter(ip, 'EpochSpec', struct('unit','min','start',2,'end',27), @(x) isempty(x) || isnumeric(x) || isstruct(x));

addParameter(ip, 'FiniteFracThresh', 0.95, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'UseZScore', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'BandHz', [0.01 0.10], @(x) isnumeric(x) && numel(x)==2);
addParameter(ip, 'WelchWinSec', 120, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'WelchOverlapFrac', 0.5, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'LDSReg', 1e-3, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'MinEigMag', 0.90, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'MaxEigMag', 1.05, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'PhaseAmpQuantileMin', 0.20, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'PLVThreshold', 0.70, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'WinSec', 600, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'StepSec', 120, @(x) isnumeric(x) && isscalar(x));

addParameter(ip, 'EnableAcfPeakFilter', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'AcfPeakLagRangeSec', [10 120], @(x) isnumeric(x) && numel(x)==2);
addParameter(ip, 'AcfMinProminence', 0.08, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'MinNumAcPeaks', 2, @(x) isnumeric(x) && isscalar(x));
addParameter(ip, 'AcfMaxLagSec', 180, @(x) isnumeric(x) && isscalar(x));

addParameter(ip, 'MakePlots', true, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'SavePlots', false, @(x) islogical(x) || isnumeric(x));
addParameter(ip, 'PlotOutDir', '', @(x) ischar(x) || isstring(x));

addParameter(ip, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

parse(ip, baseFolder, varargin{:});
P = ip.Results;

numPCsList = unique(P.NumPCsList(:)');
numPCsList = sort(numPCsList);

%% ====================== FIND VALID RECORDINGS ===========================
allFolders   = strsplit(genpath(char(baseFolder)), pathsep);
validFolders = {};
resultsFiles = {};
metaFiles    = {};
recordingID  = {};

for i = 1:numel(allFolders)
    if isempty(allFolders{i}), continue; end
    rf = fullfile(allFolders{i}, char(P.ResultsFileName));
    mf = fullfile(allFolders{i}, char(P.MetaFileName));
    if exist(rf,'file')==2 && exist(mf,'file')==2
        validFolders{end+1} = allFolders{i}; %#ok<SAGROW>
        resultsFiles{end+1} = rf; %#ok<SAGROW>
        metaFiles{end+1}    = mf; %#ok<SAGROW>
        recordingID{end+1}  = allFolders{i}; %#ok<SAGROW>
    end
end

nRec = numel(validFolders);
if nRec == 0
    error('No valid folders found under baseFolder.');
end

if P.Verbose
    fprintf('Found %d valid recordings.\n', nRec);
end

%% ====================== PRELOAD RECORDINGS =============================
REC = struct([]);
for iRec = 1:nRec
    R = load(resultsFiles{iRec});
    [dff, fps] = local_extract_dff_and_fps(R, P.DffVarCandidates, P.FpsVarCandidates, resultsFiles{iRec});
    stimFrames = local_extract_stim_frames(R, P.StimVarCandidates);
    stimFrames = unique(round(stimFrames(:)));
    stimFrames = stimFrames(isfinite(stimFrames));
    stimFrames = stimFrames(stimFrames >= 1 & stimFrames <= size(dff,2));

    epochFrames = local_resolve_epoch_spec(P.EpochSpec, fps, size(dff,2));

    REC(iRec).id         = recordingID{iRec};
    REC(iRec).folder     = validFolders{iRec};
    REC(iRec).resultsFile= resultsFiles{iRec};
    REC(iRec).metaFile   = metaFiles{iRec};
    REC(iRec).dff        = dff;
    REC(iRec).fps        = fps;
    REC(iRec).stimFrames = stimFrames;
    REC(iRec).epochFrames= epochFrames;
end

%% ====================== RUN DIMENSIONALITY SWEEP =======================
ALL = struct();
ALL.perK = struct([]);

prevAssignments = cell(nRec,1);
prevKeepIdx      = cell(nRec,1);
prevK            = NaN;

for iK = 1:numel(numPCsList)
    K = numPCsList(iK);

    if P.Verbose
        fprintf('\n====================================================\n');
        fprintf('Running numPCs = %d\n', K);
    end

    recRes = struct([]);
    rowCell = cell(nRec,1);

    for iRec = 1:nRec
        if P.Verbose
            fprintf('  Recording %d / %d\n', iRec, nRec);
        end

        OUTi = oscillator_test_one_recording_v2( ...
            REC(iRec).dff, REC(iRec).fps, ...
            'EpochFrames', REC(iRec).epochFrames, ...
            'EventFrames', REC(iRec).stimFrames, ...
            'FiniteFracThresh', P.FiniteFracThresh, ...
            'UseZScore', P.UseZScore, ...
            'NumPCs', K, ...
            'BandHz', P.BandHz, ...
            'WelchWinSec', P.WelchWinSec, ...
            'WelchOverlapFrac', P.WelchOverlapFrac, ...
            'LDSReg', P.LDSReg, ...
            'MinEigMag', P.MinEigMag, ...
            'MaxEigMag', P.MaxEigMag, ...
            'PhaseAmpQuantileMin', P.PhaseAmpQuantileMin, ...
            'PLVThreshold', P.PLVThreshold, ...
            'WinSec', P.WinSec, ...
            'StepSec', P.StepSec, ...
            'EnableAcfPeakFilter', P.EnableAcfPeakFilter, ...
            'AcfPeakLagRangeSec', P.AcfPeakLagRangeSec, ...
            'AcfMinProminence', P.AcfMinProminence, ...
            'MinNumAcPeaks', P.MinNumAcPeaks, ...
            'AcfMaxLagSec', P.AcfMaxLagSec, ...
            'MakeModeMembershipPlot', false, ...
            'DoPlots', false, ...
            'Verbose', false);

        % ---------- per-recording metrics ----------
        expRetained = sum(OUTi.pca.explained, 'omitnan');

        occModes = nnz(OUTi.mode.modeCount > 0);

        % geometric separation: average pairwise centroid distance /
        % within-mode spread in PC1-PC2-PC3
        geom = local_geometry_separation(OUTi);

        % assignment stability vs previous K
        if iK == 1
            assignStability = NaN;
            keepOverlapFrac = NaN;
        else
            [assignStability, keepOverlapFrac] = local_assignment_stability( ...
                prevKeepIdx{iRec}, prevAssignments{iRec}, ...
                OUTi.data.keepIdx, OUTi.mode.modeID);
        end

        row = struct();
        row.recordingLabel          = string(REC(iRec).id);
        row.numPCs                  = K;
        row.nKeptNeurons            = OUTi.summary.nKeptNeurons;
        row.expVarRetained          = expRetained;
        row.occupiedModes           = occModes;
        row.meanPLV                 = OUTi.summary.meanPLV;
        row.nPhaseComponents        = OUTi.summary.nPhaseComponents;
        row.nOscillatorPairsLDS     = OUTi.summary.nOscillatorPairsLDS;
        row.meanSlidingPairs        = OUTi.summary.meanSlidingPairs;
        row.maxSlidingPairs         = OUTi.summary.maxSlidingPairs;
        row.pcFreqCV                = OUTi.summary.pcFreqCV;
        row.geomCentroidSep         = geom.centroidSepMean;
        row.geomWithinSpread        = geom.withinSpreadMean;
        row.geomSepRatio            = geom.sepRatio;
        row.assignStabilityPrevK    = assignStability;
        row.keepOverlapFracPrevK    = keepOverlapFrac;

        rowCell{iRec} = row;

        recRes(iRec).OUT = OUTi; %#ok<AGROW>

        % store for next-K comparison
        prevAssignments{iRec} = OUTi.mode.modeID;
        prevKeepIdx{iRec}     = OUTi.data.keepIdx;
    end

    T = struct2table([rowCell{:}]);

    % ---------- group summary ----------
    G = struct();
    G.numPCs = K;
    G.nRecordings = height(T);

    G.meanExpVarRetained   = mean(T.expVarRetained, 'omitnan');
    G.medianExpVarRetained = median(T.expVarRetained, 'omitnan');

    G.meanOccupiedModes    = mean(T.occupiedModes, 'omitnan');
    G.medianOccupiedModes  = median(T.occupiedModes, 'omitnan');
    G.fracAllModesOccupied = mean(T.occupiedModes == K, 'omitnan');

    G.meanPLV              = mean(T.meanPLV, 'omitnan');
    G.medianPLV            = median(T.meanPLV, 'omitnan');

    G.meanNPhaseComponents = mean(T.nPhaseComponents, 'omitnan');
    G.medianNPhaseComponents = median(T.nPhaseComponents, 'omitnan');

    G.fracWithLDSOscPair   = mean(T.nOscillatorPairsLDS > 0, 'omitnan');
    G.meanLDSOscPairs      = mean(T.nOscillatorPairsLDS, 'omitnan');

    G.fracWithTransientOsc = mean(T.maxSlidingPairs > 0, 'omitnan');
    G.meanSlidingPairs     = mean(T.meanSlidingPairs, 'omitnan');

    G.meanPcFreqCV         = mean(T.pcFreqCV, 'omitnan');
    G.medianPcFreqCV       = median(T.pcFreqCV, 'omitnan');

    G.meanGeomSepRatio     = mean(T.geomSepRatio, 'omitnan');
    G.medianGeomSepRatio   = median(T.geomSepRatio, 'omitnan');

    G.meanAssignmentStabilityPrevK = mean(T.assignStabilityPrevK, 'omitnan');
    G.medianAssignmentStabilityPrevK = median(T.assignStabilityPrevK, 'omitnan');
    G.meanKeepOverlapPrevK = mean(T.keepOverlapFracPrevK, 'omitnan');

    % heuristic interpretation
    G.note = local_interpret_group_summary(K, G, prevK);

    ALL.perK(iK).K = K;
    ALL.perK(iK).recordingTable = T;
    ALL.perK(iK).groupSummary   = G;
    ALL.perK(iK).recordings     = recRes;

    prevK = K;
end

%% ====================== BUILD COMPARISON TABLES ========================
groupRows = [];
for iK = 1:numel(ALL.perK)
    G = ALL.perK(iK).groupSummary;
    row = struct();
    row.numPCs                    = G.numPCs;
    row.nRecordings               = G.nRecordings;
    row.meanExpVarRetained        = G.meanExpVarRetained;
    row.medianExpVarRetained      = G.medianExpVarRetained;
    row.meanOccupiedModes         = G.meanOccupiedModes;
    row.fracAllModesOccupied      = G.fracAllModesOccupied;
    row.meanPLV                   = G.meanPLV;
    row.meanNPhaseComponents      = G.meanNPhaseComponents;
    row.fracWithLDSOscPair        = G.fracWithLDSOscPair;
    row.meanLDSOscPairs           = G.meanLDSOscPairs;
    row.fracWithTransientOsc      = G.fracWithTransientOsc;
    row.meanSlidingPairs          = G.meanSlidingPairs;
    row.meanPcFreqCV              = G.meanPcFreqCV;
    row.meanGeomSepRatio          = G.meanGeomSepRatio;
    row.meanAssignmentStabilityPrevK = G.meanAssignmentStabilityPrevK;
    row.meanKeepOverlapPrevK      = G.meanKeepOverlapPrevK;
    row.note                      = string(G.note);
    groupRows = [groupRows; row]; %#ok<AGROW>
end

groupSummaryTable = struct2table(groupRows);

%% ====================== MAKE ONE DECISION TABLE ========================
decisionTable = groupSummaryTable;


% optional simple ranking
% Higher is better for:
%   expVarRetained, meanPLV, fracWithLDSOscPair, fracWithTransientOsc, geomSepRatio, assignmentStability
% Lower is better for:
%   meanNPhaseComponents, meanPcFreqCV
decisionTable.score = nan(height(decisionTable),1);

if height(decisionTable) >= 2
    zExp   = local_z(decisionTable.meanExpVarRetained);
    zPLV   = local_z(decisionTable.meanPLV);
    zLDS   = local_z(decisionTable.fracWithLDSOscPair);
    zTr    = local_z(decisionTable.fracWithTransientOsc);
    zGeom  = local_z(decisionTable.meanGeomSepRatio);
    zStab  = local_z(decisionTable.meanAssignmentStabilityPrevK);

    zPhaseComp = local_z(-decisionTable.meanNPhaseComponents);
    zFreqCV    = local_z(-decisionTable.meanPcFreqCV);

    decisionTable.score = ...
        zExp + zPLV + zLDS + zTr + zGeom + zStab + zPhaseComp + zFreqCV;
end


%% ====================== MAKE GROUP-LEVEL PLOTS ========================
fig = struct();

if P.MakePlots
    fig = make_numPCs_comparison_plots(groupSummaryTable, decisionTable, ALL);

    if P.SavePlots
        plotDir = char(P.PlotOutDir);
        if isempty(plotDir)
            plotDir = fullfile(char(baseFolder), 'numPCs_group_summary_plots');
        end
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end

        save_fig_if_valid(fig, 'overview',   plotDir, 'numPCs_overview.png');
        save_fig_if_valid(fig, 'stability',  plotDir, 'numPCs_stability.png');
        save_fig_if_valid(fig, 'recordings', plotDir, 'numPCs_recording_heatmaps.png');
        save_fig_if_valid(fig, 'decision',   plotDir, 'numPCs_decision.png');
    end
end
%% ====================== PRINT COMPACT OUTPUT ===========================
if P.Verbose
    fprintf('\n================ GROUP SUMMARY BY numPCs ================\n');
    disp(groupSummaryTable);

    fprintf('\n================ DECISION TABLE =========================\n');
    disp(sortrows(decisionTable, 'score', 'descend'));
end

OUT = struct();
OUT.params = P;
OUT.validFolders = validFolders;
OUT.resultsFiles = resultsFiles;
OUT.metaFiles = metaFiles;
OUT.recordingCache = REC;
OUT.all = ALL;
OUT.groupSummaryTable = groupSummaryTable;
OUT.decisionTable = decisionTable;
OUT.fig = fig;

end

%% ======================================================================
function [dff, fps] = local_extract_dff_and_fps(R, dffCandidates, fpsCandidates, filePath)
dff = [];
fps = [];

vars = fieldnames(R);

for i = 1:numel(dffCandidates)
    vn = dffCandidates{i};
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

for i = 1:numel(fpsCandidates)
    vn = fpsCandidates{i};
    if isfield(R, vn) && isnumeric(R.(vn)) && isscalar(R.(vn))
        fps = double(R.(vn));
        break;
    end
end

if isempty(fps)
    for i = 1:numel(vars)
        x = R.(vars{i});
        if isstruct(x)
            for j = 1:numel(fpsCandidates)
                vn = fpsCandidates{j};
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

%% ======================================================================
function stimFrames = local_extract_stim_frames(S, stimCandidates)
stimFrames = [];
vars = fieldnames(S);

for i = 1:numel(stimCandidates)
    vn = stimCandidates{i};
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
        for j = 1:numel(stimCandidates)
            vn = stimCandidates{j};
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

%% ======================================================================
function epochFrames = local_resolve_epoch_spec(epochSpec, fps, nFramesTotal)
if isempty(epochSpec)
    epochFrames = [1 nFramesTotal];
    return;
end

if isnumeric(epochSpec) && numel(epochSpec) == 2
    epochFrames = round(epochSpec(:))';
    epochFrames(1) = max(1, epochFrames(1));
    epochFrames(2) = min(nFramesTotal, epochFrames(2));
    if epochFrames(2) <= epochFrames(1)
        error('Resolved epoch frame interval is invalid.');
    end
    return;
end

if isstruct(epochSpec)
    if ~isfield(epochSpec,'unit') || ~isfield(epochSpec,'start') || ~isfield(epochSpec,'end')
        error('epochSpec struct must contain unit/start/end.');
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

    epochFrames = [max(1,startFrame), min(nFramesTotal,endFrame)];
    if epochFrames(2) <= epochFrames(1)
        error('Resolved epoch frame interval is invalid.');
    end
    return;
end

error('epochSpec must be [], numeric pair, or struct.');
end

%% ======================================================================
function geom = local_geometry_separation(OUT)
geom = struct('centroidSepMean',NaN,'withinSpreadMean',NaN,'sepRatio',NaN);

latentTS = OUT.pca.latentTS;
domMode_t = OUT.mode.domMode_t(:);
K = OUT.summary.numPCs;

if size(latentTS,2) >= 3
    Y = latentTS(:,1:3);
else
    Y = latentTS(:,1:min(2,size(latentTS,2)));
end

centroids = nan(K, size(Y,2));
withinSpread = nan(K,1);

for k = 1:K
    idx = (domMode_t == k);
    if nnz(idx) < 5
        continue;
    end
    Yk = Y(idx,:);
    centroids(k,:) = mean(Yk,1,'omitnan');
    diffs = Yk - centroids(k,:);
    withinSpread(k) = mean(sqrt(sum(diffs.^2,2)), 'omitnan');
end

validC = all(isfinite(centroids),2);
C = centroids(validC,:);
if size(C,1) >= 2
    D = pdist(C);
    centroidSepMean = mean(D, 'omitnan');
else
    centroidSepMean = NaN;
end

withinSpreadMean = mean(withinSpread(isfinite(withinSpread)), 'omitnan');
sepRatio = centroidSepMean / withinSpreadMean;

geom.centroidSepMean = centroidSepMean;
geom.withinSpreadMean = withinSpreadMean;
geom.sepRatio = sepRatio;
end

%% ======================================================================
function [stability, keepOverlapFrac] = local_assignment_stability(prevKeepIdx, prevModeID, keepIdx, modeID)
% Compare assignments only on neurons kept in both runs.
[commonIdx, iaPrev, iaNow] = intersect(prevKeepIdx(:), keepIdx(:), 'stable');

if isempty(commonIdx)
    stability = NaN;
    keepOverlapFrac = 0;
    return;
end

mPrev = prevModeID(iaPrev);
mNow  = modeID(iaNow);

% raw label matching is acceptable here because modes are ordered PCs, not arbitrary kmeans labels
stability = mean(mPrev == mNow, 'omitnan');
keepOverlapFrac = numel(commonIdx) / max(numel(unique([prevKeepIdx(:); keepIdx(:)])), 1);
end

%% ======================================================================
function z = local_z(x)
x = double(x);
mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');
if ~isfinite(sd) || sd == 0
    z = zeros(size(x));
else
    z = (x - mu) ./ sd;
end
end

%% ======================================================================
function txt = local_interpret_group_summary(K, G, prevK)
if isnan(prevK)
    txt = sprintf('Baseline run for K=%d.', K);
    return;
end

bits = {};

if G.meanAssignmentStabilityPrevK > 0.85
    bits{end+1} = 'assignments stable vs previous K'; %#ok<AGROW>
elseif G.meanAssignmentStabilityPrevK > 0.65
    bits{end+1} = 'assignments moderately stable vs previous K'; %#ok<AGROW>
else
    bits{end+1} = 'assignments shift substantially vs previous K'; %#ok<AGROW>
end

if G.fracAllModesOccupied > 0.9
    bits{end+1} = 'most recordings occupy all requested modes'; %#ok<AGROW>
else
    bits{end+1} = 'some requested modes remain unoccupied'; %#ok<AGROW>
end

if G.meanGeomSepRatio > 1.5
    bits{end+1} = 'good geometric separation'; %#ok<AGROW>
elseif G.meanGeomSepRatio > 1.0
    bits{end+1} = 'moderate geometric separation'; %#ok<AGROW>
else
    bits{end+1} = 'weak geometric separation'; %#ok<AGROW>
end

if G.fracWithLDSOscPair > 0.25
    bits{end+1} = 'oscillator-pair detection relatively common'; %#ok<AGROW>
else
    bits{end+1} = 'oscillator-pair detection remains sparse'; %#ok<AGROW>
end

txt = strjoin(bits, '; ');
end

%% ======================================================================
function fig = make_numPCs_comparison_plots(groupSummaryTable, decisionTable, ALL)
% MAKE_NUMPCS_COMPARISON_PLOTS
% Create compact visual summaries of how the metrics change with K.

fig = struct();

K = groupSummaryTable.numPCs;

%% ===================== FIGURE 1: OVERVIEW CURVES ======================
fig.overview = figure('Color','w', 'Position', [100 100 1500 900], ...
    'Name', 'numPCs comparison overview');

tlo = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

% 1) assignment stability + occupied modes
ax1 = nexttile(tlo); hold(ax1, 'on');
plot(ax1, K, groupSummaryTable.meanAssignmentStabilityPrevK, '-o', 'LineWidth', 1.8);
plot(ax1, K, groupSummaryTable.fracAllModesOccupied, '-s', 'LineWidth', 1.8);
xlabel(ax1, 'numPCs');
ylabel(ax1, 'value');
title(ax1, 'Stability and occupancy');
legend(ax1, {'Assignment stability', 'Frac. all modes occupied'}, 'Location','best');
grid(ax1, 'on');

% 2) geometric separation + PC freq CV
ax2 = nexttile(tlo); hold(ax2, 'on');
plot(ax2, K, groupSummaryTable.meanGeomSepRatio, '-o', 'LineWidth', 1.8);
plot(ax2, K, groupSummaryTable.meanPcFreqCV, '-s', 'LineWidth', 1.8);
xlabel(ax2, 'numPCs');
ylabel(ax2, 'value');
title(ax2, 'Geometry and timescale spread');
legend(ax2, {'Geom. sep ratio', 'Mean PC freq CV'}, 'Location','best');
grid(ax2, 'on');

% 3) oscillator detection
ax3 = nexttile(tlo); hold(ax3, 'on');
plot(ax3, K, groupSummaryTable.fracWithLDSOscPair, '-o', 'LineWidth', 1.8);
plot(ax3, K, groupSummaryTable.fracWithTransientOsc, '-s', 'LineWidth', 1.8);
plot(ax3, K, groupSummaryTable.meanSlidingPairs, '-d', 'LineWidth', 1.8);
xlabel(ax3, 'numPCs');
ylabel(ax3, 'value');
title(ax3, 'Oscillator detection');
legend(ax3, {'Frac. with LDS pair', 'Frac. with transient osc', 'Mean sliding pairs'}, 'Location','best');
grid(ax3, 'on');

% 4) PLV + phase components
ax4 = nexttile(tlo); hold(ax4, 'on');
plot(ax4, K, groupSummaryTable.meanPLV, '-o', 'LineWidth', 1.8);
plot(ax4, K, groupSummaryTable.meanNPhaseComponents, '-s', 'LineWidth', 1.8);
xlabel(ax4, 'numPCs');
ylabel(ax4, 'value');
title(ax4, 'Phase structure');
legend(ax4, {'Mean PLV', 'Mean # phase components'}, 'Location','best');
grid(ax4, 'on');

% 5) decision score
ax5 = nexttile(tlo); hold(ax5, 'on');
plot(ax5, decisionTable.numPCs, decisionTable.score, '-o', 'LineWidth', 2);
[~, bestIdx] = max(decisionTable.score);
plot(ax5, decisionTable.numPCs(bestIdx), decisionTable.score(bestIdx), 'p', 'MarkerSize', 14, 'LineWidth', 2);
xlabel(ax5, 'numPCs');
ylabel(ax5, 'score');
title(ax5, 'Aggregate decision score');
grid(ax5, 'on');

% 6) compact text summary
ax6 = nexttile(tlo);
axis(ax6, 'off');

bestK = decisionTable.numPCs(bestIdx);
txt = {
    sprintf('Best K by score: %d', bestK)
    ' '
    'Interpretation:'
    '- Low K underfits / merges structure'
    '- Higher K improves geometry and'
    '  transient oscillator sensitivity'
    '- But too high K reduces stability'
    '  and may leave modes unoccupied'
    ' '
    'Use K where:'
    '- assignments are stable'
    '- occupancy stays high'
    '- geometry improves'
    '- oscillator sensitivity improves'
    '- without obvious fragmentation'
    };
text(ax6, 0.02, 0.98, txt, 'VerticalAlignment','top', 'FontName','Courier');

%% ===================== FIGURE 2: STABILITY-FOCUSED ====================
fig.stability = figure('Color','w', 'Position', [120 120 1200 850], ...
    'Name', 'numPCs stability and tradeoffs');

tlo2 = tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

% assignment stability
ax1 = nexttile(tlo2);
bar(ax1, K, groupSummaryTable.meanAssignmentStabilityPrevK);
xlabel(ax1, 'numPCs');
ylabel(ax1, 'mean assignment stability vs previous K');
title(ax1, 'Neuron assignment stability');

% occupancy
ax2 = nexttile(tlo2);
bar(ax2, K, groupSummaryTable.meanOccupiedModes);
hold(ax2, 'on');
plot(ax2, K, K, 'k--', 'LineWidth', 1.2);
xlabel(ax2, 'numPCs');
ylabel(ax2, 'mean occupied modes');
title(ax2, 'Occupied modes vs requested modes');
legend(ax2, {'Observed', 'Ideal = K'}, 'Location','best');

% geometry vs stability
ax3 = nexttile(tlo2); hold(ax3, 'on');
scatter(ax3, groupSummaryTable.meanGeomSepRatio, groupSummaryTable.meanAssignmentStabilityPrevK, 80, K, 'filled');
for i = 1:numel(K)
    text(groupSummaryTable.meanGeomSepRatio(i), groupSummaryTable.meanAssignmentStabilityPrevK(i), ...
        sprintf('  K=%d', K(i)), 'FontSize', 9);
end
xlabel(ax3, 'mean geometric separation ratio');
ylabel(ax3, 'mean assignment stability');
title(ax3, 'Tradeoff: separation vs stability');
grid(ax3, 'on');
colorbar(ax3);

% transient oscillator sensitivity vs stability
ax4 = nexttile(tlo2); hold(ax4, 'on');
scatter(ax4, groupSummaryTable.fracWithTransientOsc, groupSummaryTable.meanAssignmentStabilityPrevK, 80, K, 'filled');
for i = 1:numel(K)
    text(groupSummaryTable.fracWithTransientOsc(i), groupSummaryTable.meanAssignmentStabilityPrevK(i), ...
        sprintf('  K=%d', K(i)), 'FontSize', 9);
end
xlabel(ax4, 'fraction with transient oscillators');
ylabel(ax4, 'mean assignment stability');
title(ax4, 'Tradeoff: transient sensitivity vs stability');
grid(ax4, 'on');
colorbar(ax4);

%% ===================== FIGURE 3: RECORDING-LEVEL HEATMAPS ============
fig.recordings = figure('Color','w', 'Position', [140 140 1600 900], ...
    'Name', 'numPCs recording-level heatmaps');

% Build recording x K matrices
Klist = [ALL.perK.K];
nK = numel(Klist);
recNames = ALL.perK(1).recordingTable.recordingLabel;
nRec = numel(recNames);

mat_stab   = nan(nRec, nK);
mat_geom   = nan(nRec, nK);
mat_plv    = nan(nRec, nK);
mat_occ    = nan(nRec, nK);
mat_pairs  = nan(nRec, nK);
mat_slide  = nan(nRec, nK);

for iK = 1:nK
    T = ALL.perK(iK).recordingTable;
    mat_stab(:,iK)  = T.assignStabilityPrevK;
    mat_geom(:,iK)  = T.geomSepRatio;
    mat_plv(:,iK)   = T.meanPLV;
    mat_occ(:,iK)   = T.occupiedModes;
    mat_pairs(:,iK) = T.nOscillatorPairsLDS;
    mat_slide(:,iK) = T.maxSlidingPairs;
end

tlo3 = tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile(tlo3);
imagesc(ax1, Klist, 1:nRec, mat_stab);
xlabel(ax1, 'numPCs');
ylabel(ax1, 'recording');
title(ax1, 'Assignment stability');
colorbar(ax1);

ax2 = nexttile(tlo3);
imagesc(ax2, Klist, 1:nRec, mat_geom);
xlabel(ax2, 'numPCs');
ylabel(ax2, 'recording');
title(ax2, 'Geometric separation');
colorbar(ax2);

ax3 = nexttile(tlo3);
imagesc(ax3, Klist, 1:nRec, mat_plv);
xlabel(ax3, 'numPCs');
ylabel(ax3, 'recording');
title(ax3, 'Mean PLV');
colorbar(ax3);

ax4 = nexttile(tlo3);
imagesc(ax4, Klist, 1:nRec, mat_occ);
xlabel(ax4, 'numPCs');
ylabel(ax4, 'recording');
title(ax4, 'Occupied modes');
colorbar(ax4);

ax5 = nexttile(tlo3);
imagesc(ax5, Klist, 1:nRec, mat_pairs);
xlabel(ax5, 'numPCs');
ylabel(ax5, 'recording');
title(ax5, 'Whole-recording LDS pairs');
colorbar(ax5);

ax6 = nexttile(tlo3);
imagesc(ax6, Klist, 1:nRec, mat_slide);
xlabel(ax6, 'numPCs');
ylabel(ax6, 'recording');
title(ax6, 'Max sliding-window pairs');
colorbar(ax6);

%% ===================== FIGURE 4: DECISION SUMMARY =====================
fig.decision = figure('Color','w', 'Position', [160 160 1200 700], ...
    'Name', 'numPCs decision summary');

ax = axes(fig.decision);
axis(ax, 'off');

[~, ord] = sort(decisionTable.score, 'descend');

lines = cell(height(decisionTable)+6,1);
lines{1} = 'Recommended dimensionality ranking';
lines{2} = ' ';
for i = 1:height(decisionTable)
    ii = ord(i);
    lines{2+i} = sprintf('%d) K=%d | score=%.3f | %s', ...
        i, decisionTable.numPCs(ii), decisionTable.score(ii), char(decisionTable.note(ii)));
end

lines{height(decisionTable)+4} = ' ';
lines{height(decisionTable)+5} = 'Rule of thumb:';
lines{height(decisionTable)+6} = 'If K=4–6 look similar, choose 6 as the descriptive working model; if K=8 adds little but reduces stability, it is probably too large.';

text(ax, 0.02, 0.98, lines, 'VerticalAlignment','top', 'FontName','Courier');
end

%% ======================================================================
function save_fig_if_valid(figStruct, fieldName, outDir, fileName)
if isfield(figStruct, fieldName)
    h = figStruct.(fieldName);
    if ~isempty(h) && isgraphics(h)
        exportgraphics(h, fullfile(outDir, fileName), 'Resolution', 200);
    end
end
end