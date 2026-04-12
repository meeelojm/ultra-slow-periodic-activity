function [U, Te] = build_synced_uni_traces(S, baseRoots, fish, cond)
% Build synchronized & cropped ethogram traces on a uniform, imaging-aligned clock.
% Returns:
%   U  : struct with fields {fra_tim_uni, fra_rat, hea_uni, ope_uni, mou_uni, eye_uni}
%   Te : the uniform time vector (seconds)

    % ---------- resolve imaging & metadata files ----------
    condAliases = {cond, [cond 'act']};
    imaMat = find_latest(baseRoots, fish, condAliases, fullfile('suite2p','*_0_0.mat'));
    metaMat = '';
    for fn = {'metadata_multimodal.mat','metadata.mat'}
        try
            metaMat = find_latest(baseRoots, fish, condAliases, fullfile('suite2p', fn{1}));
            break;
        catch
        end
    end
    if isempty(metaMat)
        error('build_synced_uni_traces:meta', 'No metadata file found for %s/%s', fish, cond);
    end

    % ---------- obtain man_sec / ini_del ----------
    man_sec = []; ini_del = [];
    try
        M = load(metaMat, 'man_sec','ini_del');
        if isfield(M,'man_sec'), man_sec = M.man_sec; end
        if isfield(M,'ini_del'), ini_del = M.ini_del; end
    catch
    end
    if isempty(man_sec) || isempty(ini_del)
        try
            [man_sec, ini_del] = ext_syn_par(imaMat, metaMat);
        catch
            warning('build_synced_uni_traces:sync', ...
                'man_sec/ini_del unavailable; using identity (1,0).');
            man_sec = 1; ini_del = 0;
        end
    end

    % ---------- imaging frame times (synchronized) ----------
    fra_tim_raw = ext_tem_fra_dat(imaMat);          % imaging timestamps (seconds, raw)
    fra_tim_syn = fra_tim_raw * man_sec + ini_del;

    % ---------- strict window [sta_tim, end_tim] ----------
    sta_tim = []; end_tim = [];
    try
        MM = load(metaMat,'sta_tim','end_tim');
        if isfield(MM,'sta_tim'), sta_tim = MM.sta_tim; end
        if isfield(MM,'end_tim'), end_tim = MM.end_tim; end
    catch
    end
    if (isempty(sta_tim) || isempty(end_tim)) && isstruct(S)
        if isfield(S,'sta_tim'), sta_tim = S.sta_tim; end
        if isfield(S,'end_tim'), end_tim = S.end_tim; end
    end
    if isempty(sta_tim) || isempty(end_tim)
        warning('build_synced_uni_traces:sta_end', ...
            'sta_tim/end_tim missing; using imaging span as crop window.');
        sta_tim = fra_tim_syn(1);
        end_tim = fra_tim_syn(end);
    end

    % ---------- derive imaging-aligned frame rate over strict window ----------
    in_strict = (fra_tim_syn > sta_tim) & (fra_tim_syn <= end_tim);
    if nnz(in_strict) >= 2
        dt = median(diff(fra_tim_syn(in_strict)));
        fra_rat = 1/max(dt, eps);
    else
        dt = median(diff(fra_tim_syn));
        fra_rat = 1/max(dt, eps);
    end

    % ---------- helper: pick best available series + its timebase ----------
    function [t_src, y_src] = pick_series(S, base) %#ok<INUSL>
        t_src = []; y_src = [];
        sigField1 = [base '_sig'];         % e.g., 'hea_sig'
        sigField2 = [base '_uni_frames'];  % e.g., 'hea_uni_frames'

        if isfield(S, sigField1) && ~isempty(S.(sigField1))
            y_src = S.(sigField1)(:);
            if isfield(S,'fra_tim') && numel(S.fra_tim)==numel(y_src)
                t_src = S.fra_tim(:);
                return;
            end
        end
        if isfield(S, sigField2) && ~isempty(S.(sigField2))
            y_src = S.(sigField2)(:);
            if isfield(S,'fra_tim_uni_frames') && numel(S.fra_tim_uni_frames)==numel(y_src)
                t_src = S.fra_tim_uni_frames(:);
                return;
            elseif isfield(S,'fra_tim') && numel(S.fra_tim)==numel(y_src)
                t_src = S.fra_tim(:);
                return;
            end
        end
        % if nothing matched, leave empty; caller will skip
    end

    % ---------- resample any present channel to the uniform, imaging-aligned grid ----------
    % Uniform grid on [sta_tim,end_tim] with imaging-derived rate:
    Te = (sta_tim + (0:floor((end_tim-sta_tim)*fra_rat))'/fra_rat);

    resamp_uni = @(t,y) smo_uni2(t, y, sta_tim, end_tim, fra_rat);

    % HEART
    [t_h, y_h] = pick_series(S, 'hea');
    if ~isempty(t_h)
        [hea_uni, Te] = resamp_uni(t_h, y_h);
    else
        hea_uni = [];
    end

    % OPERCULUM
    [t_o, y_o] = pick_series(S, 'ope');
    if ~isempty(t_o)
        [ope_uni, ~] = resamp_uni(t_o, y_o);
    else
        ope_uni = [];
    end

    % MOUTH
    [t_m, y_m] = pick_series(S, 'mou');
    if ~isempty(t_m)
        [mou_uni, ~] = resamp_uni(t_m, y_m);
    else
        mou_uni = [];
    end

    % EYE
    [t_e, y_e] = pick_series(S, 'eye');
    if ~isempty(t_e)
        [eye_uni, ~] = resamp_uni(t_e, y_e);
    else
        eye_uni = [];
    end

    % ---------- pack ----------
    U = struct();
    U.fra_tim_uni = Te(:);
    U.fra_rat     = fra_rat;
    U.hea_uni     = hea_uni(:);
    U.ope_uni     = ope_uni(:);
    U.mou_uni     = mou_uni(:);
    U.eye_uni     = eye_uni(:);
end
