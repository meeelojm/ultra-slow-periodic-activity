function RES = analyze_single_recording(inFile, outFolder, cfg)
% ANALYZE_SINGLE_RECORDING
%
% Loads one recording and performs the early preprocessing steps:
%   1) load / align variables
%   2) z-score neurons
%   3) define quiet vs active windows from bout rate
%   4) classify oscillatory neurons
%
% Later sections will append the advanced analyses and figures.

% -------------------------------------------------------------------------
% LOAD
% -------------------------------------------------------------------------
S = load(inFile, 'fps', 'dff_new', 'position', 'tai_ang_uni_frames', 'bou_ons_fra');

fps   = double(S.fps(:));
fps   = fps(1);

dff   = double(S.dff_new);              % neurons x time
pos   = double(S.position);
tail  = double(S.tai_ang_uni_frames(:)');
bouts = double(S.bou_ons_fra(:)') - 1;  % convert to 0-based internal indexing
bouts = round(bouts);

[nNeurons, T] = size(dff);
recordingMin = T / fps / 60;

% -------------------------------------------------------------------------
% BASIC SANITY / ALIGNMENT
% -------------------------------------------------------------------------
% Keep only valid bout indices
bouts = bouts(bouts >= 0 & bouts < T);

% Make sure positions are neurons x 3
if size(pos,1) >= nNeurons
    posXYZ = pos(1:nNeurons, 1:min(3,size(pos,2)));
else
    posXYZ = nan(nNeurons, 3);
    posXYZ(1:size(pos,1), 1:min(3,size(pos,2))) = pos(:, 1:min(3,size(pos,2)));
end

if size(posXYZ,2) < 3
    posXYZ(:, end+1:3) = nan;
end

% Make sure tail length matches recording length as closely as possible
if numel(tail) < T
    tail(end+1:T) = nan;
elseif numel(tail) > T
    tail = tail(1:T);
end

% -------------------------------------------------------------------------
% REMOVE ONLY FULLY INVALID NEURONS
% -------------------------------------------------------------------------
% Keep neurons that are not entirely NaN
validRows = ~all(isnan(dff), 2);
dff = dff(validRows, :);
posXYZ = posXYZ(validRows, :);
nNeurons = size(dff, 1);

% -------------------------------------------------------------------------
% Z-SCORE EACH NEURON
% -------------------------------------------------------------------------
mu = nanmean(dff, 2);
sd = nanstd(dff, 0, 2);
sd(sd == 0) = 1;

Z = (dff - mu) ./ sd;   % neurons x time
net = nanmean(Z, 1);    % population mean

% -------------------------------------------------------------------------
% BUILD BOUT TRAIN
% -------------------------------------------------------------------------
boutTrain = zeros(1, T);
boutTrain(bouts + 1) = 1;   % back to MATLAB indexing

% -------------------------------------------------------------------------
% DEFINE QUIET / ACTIVE WINDOWS
% -------------------------------------------------------------------------
win  = round(cfg.winMin  * 60 * fps);
step = round(cfg.stepMin * 60 * fps);

starts = 1:step:(T - win + 1);

winTable = table();
winTable.WindowID = (1:numel(starts))';
winTable.StartFrame = starts(:);
winTable.EndFrame   = (starts(:) + win - 1);
winTable.StartMin   = (winTable.StartFrame - 1) / fps / 60;
winTable.EndMin     = (winTable.EndFrame   - 1) / fps / 60;

winTable.BoutRatePerMin = nan(height(winTable),1);
winTable.TailRMS        = nan(height(winTable),1);
winTable.TailSTD        = nan(height(winTable),1);

for i = 1:height(winTable)
    s = winTable.StartFrame(i);
    e = winTable.EndFrame(i);

    winTable.BoutRatePerMin(i) = sum(boutTrain(s:e)) / (win / fps / 60);

    tailSeg = tail(s:e);
    winTable.TailRMS(i) = rms(tailSeg, 'omitnan');
    winTable.TailSTD(i) = std(tailSeg, 0, 'omitnan');
end

[quietBoutRate, iq]  = min(winTable.BoutRatePerMin);
[activeBoutRate, ia] = max(winTable.BoutRatePerMin);

qs = winTable.StartFrame(iq);
qe = winTable.EndFrame(iq);

as = winTable.StartFrame(ia);
ae = winTable.EndFrame(ia);

quietZ  = Z(:, qs:qe);
activeZ = Z(:, as:ae);

quietTail  = tail(qs:qe);
activeTail = tail(as:ae);

% -------------------------------------------------------------------------
% WHOLE-RECORDING DOMINANT FREQUENCY / BANDPOWER
% -------------------------------------------------------------------------
domfAll = nan(nNeurons,1);
bpAll   = nan(nNeurons,1);

for n = 1:nNeurons
    [domfAll(n), bpAll(n)] = domf_bp(Z(n,:), fps, cfg.fullBand);
end

% -------------------------------------------------------------------------
% DEFINE OSCILLATORY NEURONS
% -------------------------------------------------------------------------
oscMask = domfAll >= cfg.coreBand(1) & domfAll <= cfg.coreBand(2);
oscIdx  = find(oscMask);
nOsc    = numel(oscIdx);

if nOsc >= 1
    Zosc = Z(oscIdx, :);
    popOsc = mean(Zosc, 1, 'omitnan');
else
    Zosc = [];
    popOsc = net;
end

% -------------------------------------------------------------------------
% STORE EARLY RESULTS IN RES
% -------------------------------------------------------------------------
RES = struct();

RES.file         = inFile;
RES.outFolder    = outFolder;

RES.fps          = fps;
RES.nNeurons     = nNeurons;
RES.nOsc         = nOsc;
RES.recordingMin = recordingMin;

RES.dff          = dff;
RES.Z            = Z;
RES.Zosc         = Zosc;
RES.net          = net;
RES.popOsc       = popOsc;

RES.posXYZ       = posXYZ;
RES.tail         = tail;
RES.bouts        = bouts;
RES.boutTrain    = boutTrain;

RES.winTable         = winTable;
RES.quietBoutRate    = quietBoutRate;
RES.activeBoutRate   = activeBoutRate;

RES.qs = qs;
RES.qe = qe;
RES.as = as;
RES.ae = ae;

RES.quietZ     = quietZ;
RES.activeZ    = activeZ;
RES.quietTail  = quietTail;
RES.activeTail = activeTail;

RES.domfAll = domfAll;
RES.bpAll   = bpAll;
RES.oscMask = oscMask;
RES.oscIdx  = oscIdx;

% -------------------------------------------------------------------------
% QUICK QC PRINTS
% -------------------------------------------------------------------------
fprintf('  fps: %.3f Hz\n', fps);
fprintf('  neurons: %d\n', nNeurons);
fprintf('  recording length: %.2f min\n', recordingMin);
fprintf('  oscillatory neurons (%.3f-%.3f Hz): %d\n', cfg.coreBand(1), cfg.coreBand(2), nOsc);
fprintf('  quiet window: %.2f-%.2f min, %.3f bouts/min\n', ...
    (qs-1)/fps/60, (qe-1)/fps/60, quietBoutRate);
fprintf('  active window: %.2f-%.2f min, %.3f bouts/min\n', ...
    (as-1)/fps/60, (ae-1)/fps/60, activeBoutRate);

% -------------------------------------------------------------------------
% PLACEHOLDER FOR NEXT SECTIONS
% -------------------------------------------------------------------------
% Section 3 onward will append more fields into RES and export figures.

end
