%% ETHOGRAM -> NEURON TIME  + dominant-frequency per fish
% - For each fish/condition:
%     * load ethogram + Suite2p
%     * map ethogram signals to neuron frameTimes
%     * save ethogram_in_neurontime.mat (one file per condition)
% - Then, for each fish:
%     * load & concatenate mapped signals across conditions (in given order)
%     * run calcsort_autocorr_freq_analysis_v3 per channel (SLOW bands)
%
% You need calcsort_autocorr_freq_analysis_v3 on the path.

clear; clc; set(0,'defaultfigurecolor',[1 1 1]);

%% ---------- CONFIG ----------
baseRoots = {'\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low', ...
             'Z:\temp\Emilian temp\ccn2a\high2low', ...
             'Z:\Emilian temp\ccn2a\high2low'};


fishList  = {'f010','f011','f012'};
condList  = {'high','mid','low'};     % order = concatenation order
% Accept both "...act" and plain condition, prefer the *act variant first
condAliases = @(c){[c 'act'], c};


% calcsort parameters (for slow-mod analysis in PASS-2)
acfg.NumLags          = 1500;
acfg.windowSize       = 1024;
acfg.maxpx            = 15;
acfg.minPeakThreshold = 1;

tailSmoothSec = 0.10;  % for |dθ/dt| pre-derivative smoothing
fprintf('== Dry-run: resolving input files ==\n');
for f = 1:numel(fishList)
    fish = fishList{f};
    for c = 1:numel(condList)
        cond = condList{c};
        fprintf('Fish %s | Cond %s\n', fish, cond);
        [~, triedE] = find_latest(baseRoots, fish, condAliases(cond), fullfile('suite2p','ethogram_*.mat'), true);
        [~, triedS] = find_latest(baseRoots, fish, condAliases(cond), fullfile('suite2p','dffs_repact_respcells.mat'), true);
        % (The function already prints the chosen file; the variables are unused here.)
    end
end

%% ---------- PASS 1 (simple): map each (fish,cond) ethogram to neuron time ----------
fprintf('== Mapping ethograms to neuron-time (simple, no extra sync) ==\n');

for f = 1:numel(fishList)
  fish = fishList{f};
  for c = 1:numel(condList)
    cond = condList{c};

    ethoFile = find_latest(baseRoots, fish, condAliases(cond), fullfile('suite2p','ethogram_*.mat'));
    s2pFile  = find_latest(baseRoots, fish, condAliases(cond), fullfile('suite2p','dffs_repact_respcells.mat'));

    S = load(ethoFile);
    R = load(s2pFile, 'dff_new','fps');
    fpsN = double(R.fps);
    Tn   = (0:size(R.dff_new,2)-1)'/fpsN;        % neuron-time (s)

    % Helper: resample with safe length handling
    resamp = @(t,y) interp1(t(:), y(:), Tn, 'linear', 'extrap');

    OUT = struct();
    OUT.neuron_time_sec = Tn;
    OUT.fps_neuro       = fpsN;
    OUT.source_files    = struct('ethogram', ethoFile, 'suite2p', s2pFile);

    % -------- tail angle + tail speed (use fra_tim_uni_frames if present) --------
    if isfield(S,'tai_ang_uni_frames') && ~isempty(S.tai_ang_uni_frames)
        tA = get_timevec(S,'fra_tim_uni_frames', numel(S.tai_ang_uni_frames));
        yA = S.tai_ang_uni_frames(:);
        n  = min(numel(tA), numel(yA));  tA=tA(1:n); yA=yA(1:n);

        OUT.tail_angle = resamp(tA, yA);

        % tail_speed = |dθ/dt| computed on ethogram time, then resampled
        yAsm = movmean(yA, max(1, round(tailSmoothSec/(median(diff(tA))+eps))));
        dAdt = abs(gradient(yAsm, tA));
        OUT.tail_speed = resamp(tA, dAdt);
    end

    % -------- VR speed (uni_frames) --------
    if isfield(S,'vir_spe_uni_frames') && ~isempty(S.vir_spe_uni_frames)
        tV = get_timevec(S,'fra_tim_uni_frames', numel(S.vir_spe_uni_frames));
        yV = S.vir_spe_uni_frames(:);
        n  = min(numel(tV), numel(yV));  tV=tV(1:n); yV=yV(1:n);
        OUT.vr_speed = resamp(tV, yV);
    end

    % -------- Heart / Operculum / Mouth / Eye (on fra_tim) --------
    [tH,yH] = get_sig_from_fra(S,'hea_sig');  if ~isempty(tH), OUT.heart     = resamp(tH,yH); end
    [tO,yO] = get_sig_from_fra(S,'ope_sig');  if ~isempty(tO), OUT.operculum = resamp(tO,yO); end
    [tM,yM] = get_sig_from_fra(S,'mou_sig');  if ~isempty(tM), OUT.mouth     = resamp(tM,yM); end
    [tE,yE] = get_sig_from_fra(S,'eye_sig');  if ~isempty(tE), OUT.eye       = resamp(tE,yE); end


    % save per-condition mapping
    outFile = fullfile(fileparts(ethoFile), 'ethogram_in_neurontime.mat');
    save(outFile, '-struct', 'OUT', '-v7.3');

    % quick sanity print
    fprintf('Saved: %s | fields: %s\n', outFile, strjoin(setdiff(fieldnames(OUT), {'source_files'}), ', '));

  end
end
%% ---------- PASS 2: concatenate per fish and run calcsort on slow modulations ----------
fprintf('\n== Concatenating per fish and running calcsort ==\n');
DOM = struct();

for f = 1:numel(fishList)
  fish = fishList{f};

  % Initialize accumulator
  cat = struct('neuron_time_sec',[], ...
               'tail_angle',[], 'tail_speed',[], 'vr_speed',[], ...
               'heart',[], 'operculum',[], 'mouth',[], 'eye',[], ...
               'heart_env',[], 'operc_env',[], 'mouth_env',[], 'eye_env',[], ...
               'tail_IF_window',[], 'heart_IF_window',[], 'operc_IF_window',[], 'mouth_IF_window',[], 'eye_IF_window',[], ...
               'tail_bout_rate',[], 'heart_bout_rate',[], 'operc_bout_rate',[], 'mouth_bout_rate',[]);
  fps_ref = [];
  outdir_fish = '';

  % Load per-condition neuron-time files (in desired order) and concatenate
  for c = 1:numel(condList)
    cond  = condList{c};
    mfile = find_latest(baseRoots, fish, condAliases(cond), fullfile('suite2p','ethogram_in_neurontime.mat'));
    if c == 1
        outdir_fish = fileparts(mfile);  % existing folder to save into
    end
    M = load(mfile);

    % capture fps once
    if isempty(fps_ref) && isfield(M,'fps_neuro') && ~isempty(M.fps_neuro)
        fps_ref = double(M.fps_neuro);
    end

    % append time (make continuous)
    if isempty(cat.neuron_time_sec)
        t0 = 0;
    else
        fps_tmp = fps_ref;
        if (isempty(fps_tmp) || ~isfinite(fps_tmp) || fps_tmp<=0) && numel(cat.neuron_time_sec)>1
            fps_tmp = 1/median(diff(cat.neuron_time_sec));
        end
        if isempty(fps_tmp) || ~isfinite(fps_tmp) || fps_tmp<=0, fps_tmp = 1; end
        t0 = cat.neuron_time_sec(end) + 1/fps_tmp;
    end
    cat.neuron_time_sec = [cat.neuron_time_sec; M.neuron_time_sec(:) + t0];

    % append all known fields if present
    fns = fieldnames(cat); fns = setdiff(fns, {'neuron_time_sec'});
    for kf = 1:numel(fns)
        fn = fns{kf};
        if isfield(M, fn) && ~isempty(M.(fn))
            cat.(fn) = [cat.(fn); M.(fn)(:)];
        end
    end
  end

  % Fallback: infer fps from concatenated time if still empty
  if isempty(fps_ref) || ~isfinite(fps_ref) || fps_ref<=0
      fps_ref = 1/median(diff(cat.neuron_time_sec));
  end

  % Nyquist after fps_ref is known
  nyq = fps_ref/2;

  % SLOW modulation DF (NOT the beats)
  chNames = {'tail_speed','vr_speed','heart','operculum','mouth','eye', ...
             'tail_IF_window','heart_IF_window','operc_IF_window','mouth_IF_window','eye_IF_window', ...
             'heart_env','operc_env','mouth_env','eye_env', ...
             'tail_bout_rate','heart_bout_rate','operc_bout_rate','mouth_bout_rate'};

  DOM.(fish) = struct();
  slowBandDesired = [0 5];  % Hz

  for k = 1:numel(chNames)
    ch = chNames{k};
    if ~isfield(cat, ch) || isempty(cat.(ch)), continue; end

    x = double(cat.(ch)(:)).';

    % Nyquist clamp
    band = [max(slowBandDesired(1),0), min(slowBandDesired(2), 0.9*nyq)];
    if band(2) <= band(1) + eps
        fprintf('%s | %s: SKIP slow DF (Nyquist=%.3f Hz)\n', fish, ch, nyq);
        continue;
    end

    [acs, domf, Pxxs, Fss, isoinds] = calcsort_autocorr_freq_analysis_v3( ...
        x, acfg.NumLags, fps_ref, acfg.windowSize, acfg.maxpx, acfg.minPeakThreshold, band);

    DOM.(fish).(ch) = struct('domf',domf,'Fss',Fss,'Pxxs',Pxxs,'acs',acs,'isoinds',isoinds, ...
                             'freqRange',band,'fps',fps_ref);
    fprintf('%s | %s: slow DF = %.3f Hz (band [%g %g])\n', fish, ch, domf, band);
  end

  % Save per-fish outputs beside the mapped files
  outf = fullfile(outdir_fish, [fish '_etho_neurontime_concat_and_domfreq.mat']);
  save(outf, 'cat','DOM','-v7.3');
  fprintf('Saved per-fish concat & domfreq: %s\n', outf);
end

%% ---------------- local helpers ----------------


function [t, y] = pick_series(S, uni_field, sig_field)
% Prefer *_uni_frames if present; else *_sig with a matching time vector.
    if isfield(S,uni_field) && ~isempty(S.(uni_field))
        y = S.(uni_field)(:);
        t = pick_time_for_field(S, numel(y));
    elseif isfield(S,sig_field) && ~isempty(S.(sig_field))
        y = S.(sig_field)(:);
        t = pick_time_for_field(S, numel(y));
    else
        error('Neither %s nor %s present.', uni_field, sig_field);
    end
end

function t = pick_time_for_field(S, ylen)
% Return a time vector that matches this series length or build from fps
    cand = {'fra_tim_uni_frames','tim_uni_frames','fra_tim_uni','fra_tim','time','t'};
    t = [];
    for k = 1:numel(cand)
        nm = cand{k};
        if isfield(S, nm) && numel(S.(nm)) == ylen
            t = S.(nm)(:);
            break;
        end
    end
    if isempty(t)
        fs = [];
        if isfield(S,'fra_rat_get') && ~isempty(S.fra_rat_get), fs = double(S.fra_rat_get); end
        if isempty(fs) && isfield(S,'fra_rat') && ~isempty(S.fra_rat), fs = double(S.fra_rat); end
        assert(~isempty(fs) && fs>0, 'No matching time vector or fps found.');
        t = (0:ylen-1)'/fs;
    end
end

function yq = safe_interp_to(t, y, t_q)
% Deduplicate t and interpolate y -> t_q
    [tu, ia] = unique(t(:), 'stable');
    yu = y(ia);
    yq = interp1(tu, yu, t_q, 'linear', 'extrap');
end

function bouts = detect_bouts(t, env, fs, thr_k, min_bout_sec)
% Boolean supra-threshold with MAD-robust threshold; return [start_idx end_idx] rows
    thr = median(env,'omitnan') + thr_k*mad(env,1);
    a   = env > thr;
    a   = bwareaopen(a, max(1,round(min_bout_sec*fs)));
    L   = bwlabel(a);
    nb  = max(L);
    bouts = zeros(nb,2);
    for i=1:nb
        idx = find(L==i);
        bouts(i,:) = [idx(1), idx(end)];
    end
end

function f_list = per_bout_peakfreq(t, y, fs, bouts, band)
% Welch PSD per bout; return list of peak freqs (Hz)
    f_list = nan(0,1);
    if isempty(bouts), return; end
    for i = 1:size(bouts,1)
        seg = y(bouts(i,1):bouts(i,2));
        N   = numel(seg);
        if N < round(0.3*fs), continue; end
        nfft = 2^nextpow2(max(N,512));
        [Pxx,F] = pwelch(seg, hamming(min(N,512)), [], nfft, fs);
        m = (F>=band(1)) & (F<=band(2));
        if any(m)
            [~,ix] = max(Pxx(m));
            f = F(m); f_list(end+1,1) = f(ix); %#ok<AGROW>
        end
    end
end

function [t_win, f_win] = windowed_peakfreq(t, y, fs, win_sec, hop_sec, band)
% Sliding window Welch peak frequency; returns centers and peak Hz
    N = numel(t);
    win = max(1, round(win_sec*fs));
    hop = max(1, round(hop_sec*fs));
    if N < win, t_win = t([]); f_win = []; return; end
    idxs = 1:hop:(N-win+1);
    t_win = t(idxs + floor(win/2));
    f_win = nan(numel(idxs),1);
    nfft  = 2^nextpow2(max(win,1024));
    for k=1:numel(idxs)
        seg = y(idxs(k):idxs(k)+win-1);
        [Pxx,F] = pwelch(seg, hamming(min(numel(seg),512)), [], nfft, fs);
        m = (F>=band(1)) & (F<=band(2));
        if any(m)
            [~,ix] = max(Pxx(m));
            fw = F(m); f_win(k) = fw(ix);
        end
    end
end

function r_ts = bout_rate_ts(t, bouts, win_sec)
% Bouts/min as a continuous time series (boxcar kernel of length win_sec)
    r_ts = zeros(size(t));
    if isempty(bouts), return; end
    starts = t(bouts(:,1));
    win = win_sec;
    for i=1:numel(t)
        r_ts(i) = 60 * sum(starts >= (t(i)-win) & starts <= t(i)) / max(win,eps);
    end
end

    % ---- local helper for the 4 *_sig signals ----
    function mapSig(outName, fld)
        if isfield(S,fld) && ~isempty(S.(fld)) && isfield(S,'fra_tim') && ~isempty(S.fra_tim)
            tS = S.fra_tim(:);
            yS = S.(fld)(:);
            n  = min(numel(tS), numel(yS));  tS=tS(1:n); yS=yS(1:n);
            OUT.(outName) = resamp(tS, yS);
        end
    end


function [tS,yS] = get_sig_from_fra(S, fld)
    % Return matching time vector and signal for an ethogram *_sig field
    tS = []; yS = [];
    if isfield(S,fld) && ~isempty(S.(fld)) && isfield(S,'fra_tim') && ~isempty(S.fra_tim)
        tS = S.fra_tim(:);
        yS = S.(fld)(:);
        n  = min(numel(tS), numel(yS));
        tS = tS(1:n);
        yS = yS(1:n);
    end
end

function t = get_timevec(S, prefName, expectedLen)
    % Prefer S.(prefName) if long enough, else fall back to S.fra_tim
    if isfield(S, prefName) && numel(S.(prefName)) >= expectedLen
        t = S.(prefName)(:);
    elseif isfield(S,'fra_tim') && numel(S.fra_tim) >= expectedLen
        t = S.fra_tim(:);
    else
        error('No suitable time vector for %s (need >= %d samples).', prefName, expectedLen);
    end
end
