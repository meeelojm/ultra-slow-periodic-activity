function h = plot_ethogram_file(S, plotTitle)
%PLOT_ETHOGRAM_FILE Plot one loaded ethogram struct using plo_eth_subplots.
%
% Inputs
% ------
% S         : output of load_ethogram_file
% plotTitle : title string (optional)

    if nargin < 2 || isempty(plotTitle)
        plotTitle = sprintf('Behavior summary %s %s', S.fish, S.cond);
    end

    D = S.data;

    % ---------- required helpers ----------
    getv = @(name, default) get_field_or_default(D, name, default);

    % ---------- time limits ----------
    sta_tim = getv('sta_tim', []);
    end_tim = getv('end_tim', []);

    % ---------- event onsets ----------
    dru_ons = getv('dru_ons', []);
    sti_ons = getv('sti_ons', []);
    spo_bou_ons = getv('spo_bou_ons', []);

    % ---------- binary behavior traces ----------
    tai_bea = getv('tai_bea', []);
    hea_bea = getv('hea_bea', []);
    ope_bea = getv('ope_bea', []);
    mou_bea = getv('mou_bea', []);
    eye_bea = getv('eye_bea', []);

    % ---------- continuous signals ----------
    tai_sig = getv('tai_sig', []);
    hea_sig = getv('hea_sig', []);
    ope_sig = getv('ope_sig', []);
    mou_sig = getv('mou_sig', []);
    eye_sig = getv('eye_sig', []);

    % If you do not want to show tail continuous trace:
    tai_sig = [];

    % ---------- rate / timing ----------
    rat_tim     = getv('tim_bin', []);
    rat_tim_par = getv('rat_bin_par', []);
    fra_tim     = getv('fra_tim', []);
    tim_rob     = getv('tim_rob', []);
    rat_rob     = getv('rat_rob', []);

    % ---------- colors ----------
    col_tai = getv('col_tai', [0 0 0]);
    col_hea = getv('col_hea', [0 0 0]);
    col_ope = getv('col_ope', [0 0 0]);
    col_mou = getv('col_mou', [0 0 0]);
    col_eye = getv('col_eye', [0 0 0]);

    col_tai_sig = getv('col_tai_sig', col_tai);
    col_hea_sig = getv('col_hea_sig', col_hea);
    col_ope_sig = getv('col_ope_sig', col_ope);
    col_mou_sig = getv('col_mou_sig', col_mou);
    col_eye_sig = getv('col_eye_sig', col_eye);

    col_rob = getv('col_rob', [0 0 0]);
    col_dru = getv('col_dru', [0.5 0.5 0.5]);
    col_sti = getv('col_sti', [0.7 0.7 0.7]);
    col_spo_bou_ons = getv('col_spo_bou_ons', [0 0 0]);

    % ---------- font ----------
    fon_siz = getv('fon_siz', 10);

    % ---------- plot ----------
    h = figure('Color', 'w', ...
               'Name', sprintf('%s_%s', S.fish, S.cond), ...
               'NumberTitle', 'off');

    plo_eth_subplots( ...
        sta_tim, end_tim, ...
        dru_ons, sti_ons, tai_bea, hea_bea, ope_bea, mou_bea, eye_bea, ...
        tai_sig, hea_sig, ope_sig, mou_sig, eye_sig, ...
        rat_tim, rat_tim_par, fra_tim, tim_rob, rat_rob, spo_bou_ons, ...
        col_tai, col_hea, col_ope, col_mou, col_eye, ...
        col_tai_sig, col_hea_sig, col_ope_sig, col_mou_sig, col_eye_sig, ...
        fon_siz, col_rob, col_dru, col_sti, col_spo_bou_ons, ...
        plotTitle);

end


function v = get_field_or_default(S, fieldName, defaultValue)
% Small utility to safely read a field from a struct loaded from MAT.

    if isfield(S, fieldName)
        v = S.(fieldName);
    else
        v = defaultValue;
    end
end