%% Dominant frequencies of ethogram signals (High/Mid/Low)
% Runs your calcsort_autocorr_freq_analysis_v3 on tail/mouth/operculum/heart/VR speed

clear; clc;

% ---- roots that exist in your setup (first existing one will be used) ----
base_dirs = { ...
  '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low', ...
  'Z:\temp\Emilian temp\ccn2a\high2low', ...
  'Z:\Emilian temp\ccn2a\high2low' ...
};

fish_ids = {'f010','f011','f012'};
conds    = {'high','mid','low'};   % also accepts highact/midact/lowact

% ---- your autocorr/PSD parameters (same style you use for neurons) ----
cfg.num_lags         = 1500;
cfg.windowSize       = 1024;
cfg.maxpx            = 15;
cfg.minPeakThreshold = 1;

% ---- which signals + frequency bands (Hz) ----
siglist = [ ...
   struct('field','tai_ang_uni_frames', 'name','tail_speed',    'mode','ang', 'freqRange',[1    10 ]), ...
   struct('field','mou_sig',            'name','mouth',         'mode','raw', 'freqRange',[0.5   5 ]), ...
   struct('field','ope_sig',            'name','operculum',     'mode','raw', 'freqRange',[0.5   5 ]), ...
   struct('field','hea_sig',            'name','heart_envelope','mode','raw', 'freqRange',[2     6 ]) ...
];

smooth_sec = 0.10;  % for tail -> |dθ/dt| only

RESULTS = struct();

for s = 1:numel(siglist)
    info = siglist(s);
    fprintf('\n==== %s ====\n', info.name);
    for f = 1:numel(fish_ids)
        fish = fish_ids{f};
        for c = 1:numel(conds)
            cond = conds{c};

            fpath = find_etho_path(base_dirs, fish, cond);              % <- helper below
            [t, y, fs] = load_etho_var_robust(fpath, info.field);       % <- helper below

            % build analysis trace
            if strcmpi(info.mode,'ang')
                win  = max(1, round(smooth_sec*fs));
                y_sm = movmean(y, win, 1, 'omitnan');
                dy   = gradient(y_sm) * fs;      % dθ/dt (rad/s or a.u./s)
                trace = abs(dy);
            else
                trace = movmean(abs(y), max(1,round(0.05*fs)), 1, 'omitnan');
            end

            % call your function (expects signals x time)
            cfg.freqRange = info.freqRange;
            sig_row = double(trace(:)).';        % 1 x T
            [acs, domf, Pxxs, Fss, isoinds] = calcsort_autocorr_freq_analysis_v3( ...
                sig_row, cfg.num_lags, fs, cfg.windowSize, cfg.maxpx, cfg.minPeakThreshold, cfg.freqRange);

            RESULTS.(info.name).(fish).(cond) = struct( ...
                'domf', domf, 'Fss', Fss, 'Pxxs', Pxxs, 'acs', acs, 'isoinds', isoinds, ...
                'fs', fs, 'file', fpath);

            fprintf('%s | %s: dominant f = %.3f Hz\n', fish, cond, domf);
        end
    end
end

% ---- quick summary table per signal ----
for s = 1:numel(siglist)
    info = siglist(s);
    fprintf('\nSummary %s (Hz):\n', info.name);
    for f = 1:numel(fish_ids)
        fish = fish_ids{f};
        R = RESULTS.(info.name).(fish);
        fprintf('%s   High=%.3f   Mid=%.3f   Low=%.3f\n', fish, R.high.domf, R.mid.domf, R.low.domf);
    end
end

%% ========== Local functions (keep these after the main script) ==========



