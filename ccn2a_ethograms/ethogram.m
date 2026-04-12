%% =========================
%  0) USER SETTINGS
%  =========================
baseFolder = 'E:\Emiliano\ccn2a\high2low\20260211_09_54_55_f026_ejm_ccn2a_gcamp6s_21dpf_ong_highact\suite2p';
set(0,'defaultfigurecolor',[1 1 1])

% Tail extraction parameters
seg_len     = 1024;
mph_tai     = 0.19;
fra_mph_mpp = 0.5;
tim_bas_lim = [-1 0];
dur_tai     = 1;
lat_thr     = 0.3;
win_siz     = 120;

% Ethogram params (initial guesses; tune later)
fra_rat_get   = 120;  % Hz
pri_com_hea   = 1; pri_com_ope = 1; pri_com_mou = 1; pri_com_eye = 1;
hea_per_pro   = 70;   % was 2 -> too permissive
ope_per_pro   = 80;   % was 10
mou_per_pro   = 80;   % was 10
ope_per_wid_min = 10; ope_per_wid_max = 90; ope_per_dec = 60;
ope_dir       = 1;
mou_dir       = 0;
min_thr       = 15;   % was 1000 -> likely too high

% Plot colors (chars or RGB both OK)
col_tai='b'; col_hea='r'; col_ope='g'; col_mou='m'; col_eye='c';
col_tai_sig='b'; col_hea_sig='r'; col_ope_sig='g'; col_mou_sig='m'; col_eye_sig='c';
col_rob='k'; col_dru='k'; col_sti=[0.3 0.3 0.3]; col_spo_bou_ons='k';
fon_siz = 10;


%% =========================
%  1) FIND VALID FOLDERS
%  =========================
allFolders   = strsplit(genpath(baseFolder), pathsep);
validFolders = {};

imaPattern = '^\d{2}-\d{2}-\d{4}_\d_\d\.mat$';
robPattern = '^Robbrecht_experiment\d+_trial\d+_TID\d+\.mat$';

for i = 1:numel(allFolders)
    curF = allFolders{i};
    if isempty(curF), continue; end

    repact_file   = fullfile(curF,'dffs_repact_respcells.mat');
    metadata_file = fullfile(curF,'metadata_multimodal.mat');

    listing  = dir(fullfile(curF,'*.mat'));
    ima_list = {};
    rob_list = {};

    for k = 1:numel(listing)
        nm = listing(k).name;
        if ~isempty(regexp(nm, imaPattern,'once'))
            ima_list{end+1} = fullfile(curF, nm); %#ok<SAGROW>
        elseif ~isempty(regexp(nm, robPattern,'once'))
            rob_list{end+1} = fullfile(curF, nm); %#ok<SAGROW>
        end
    end

    if exist(repact_file,'file')==2 && exist(metadata_file,'file')==2 && ...
       ~isempty(ima_list) && ~isempty(rob_list)
        validFolders{end+1} = struct( ...
            'folder',           curF, ...
            'repact_file',      repact_file, ...
            'fil_pat_exp',      metadata_file, ...
            'fil_pat_ima_fis',  {ima_list}, ...
            'fil_pat_tai',      {rob_list} );
    end
end

fprintf('Found %d valid folders.\n', numel(validFolders));
if isempty(validFolders)
    error('No valid folders found. Check baseFolder and patterns.');
end


%% =========================
%  2) PICK A FOLDER (or loop)
%  =========================
u = 1; % or: for u = 1:numel(validFolders)

repact_file     = char(validFolders{u}.repact_file);
fil_pat_tai     = string(validFolders{u}.fil_pat_tai);      % array of paths
fil_pat_exp     = char(validFolders{u}.fil_pat_exp);
fil_pat_ima_fis = string(validFolders{u}.fil_pat_ima_fis);  % array of paths

% Load repact
S = load(repact_file, "dff_new","traces","stims","positions","fps","fishID","sti_ons","position");
if isfield(S,'position'), positions = S.position; else positions = S.positions; end
dff_new = S.dff_new; fps = S.fps; stims = S.stims; % add others if needed

% Basic time axes for this folder
nFrames    = size(dff_new,2);
frameTimes = (0:nFrames-1)'/fps;


%% =========================
%  3) EXTRACT & SYNC TAIL
%  =========================
[tai_ang, tai_ang_spe, tai_ang_uni, fra_tim_uni, vir_spe_uni, fra_tim, ...
    ~, ~, vir_spe, fra_ind, fra_rat_get_from_file, tim_tri_uni_fra, ...
    rec_tai_ang_uni_fra_sta_con, n_tri_con, n_con, tai_bea, ang_ave_con_win, ...
    nta_fra_sta_con, ave_nta_con_win, bea_lat_tri, per_pro_con, ...
    log_tri_con, bea_lat_con, ave_bas_spe, man_sec, ini_del, tai_bea_ang] ...
    = ext_tai_dat_seg(fil_pat_tai, fil_pat_exp, fil_pat_ima_fis, seg_len, ...
                      mph_tai, fra_mph_mpp, tim_bas_lim, dur_tai, lat_thr, win_siz);

% Interpolate uniform signals to imaging frame times
vir_spe_uni_frames = interp1(fra_tim_uni, vir_spe_uni, frameTimes,'linear',0);
tai_ang_uni_frames = interp1(fra_tim_uni, tai_ang_uni, frameTimes,'linear',0);

% Bout on/off in seconds and frames
[bou_ons_bou, bou_off_bou, dur_bou] = cal_bou_ons_bou(1, nFrames, tai_bea, fra_tim);
bou_ons_fra = floor(bou_ons_bou * fps);
bou_off_fra = floor(bou_off_bou * fps);

% Beat amplitudes aligned to frames (optional features)
frameIdx = round(tai_bea * fps) + 1;
frameIdx = frameIdx(frameIdx>=1 & frameIdx<=numel(frameTimes));
beatAmpImpulse = zeros(size(frameTimes));
beatAmpImpulse(frameIdx) = tai_bea_ang(1:numel(frameIdx));  % guard lengths


%% =========================
%  4) ETHOGRAM: LOAD IMAGE STACK & ROIs
%  =========================
% Choose the first imaging file as the ethogram source
imaFile = char(fil_pat_ima_fis(1));
T = load(imaFile);                % expect a variable named 'data' (HxWxT)
assert(isfield(T,'data'), 'File %s does not contain variable ''data''.', imaFile);
data = T.data;

% Crop + align full stack. For a crisp crop preview, average first 1000 frames.
fra = ext_fra(data, 1:1000);      % uses average of these for ROI preview (if your cro_ima supports it)
clear data
tifFile = fullfile(baseFolder, sprintf('registered_fra_%s.tif', datestr(now,'yyyymmdd_HHMMSS')));
save_stack_big_tiff(fra, tifFile);




%% =========================
%  5) ETHOGRAM: RUN EXTRACTION
%  =========================
% Draw ROIs on an averaged/contrast-enhanced image (as we wrote in ext_mas)
[mas_hea, mas_gil, mas_mou, mas_eye, mas_noi] = ext_mas(fra, min(10000,size(fra,3)), 'mean', true);
sta_tim = 0;
end_tim = size(fra,3) / fra_rat_get;
fra_tim = (0:size(fra,3)-1)'/fra_rat_get + sta_tim;

% No stimulus/drug schedule in this dataset
dur_spo = end_tim; n_con = 0; dur_con = 0; isi = 0; dur_con_end = 0;
tic()
[hea_bea, hea_sig, ope_bea, ope_sig, mou_bea, mou_sig, eye_bea, eye_sig, ...
 tim_bin, rat_bin_par, rat_raw_bin_par, tim_rob, rat_rob, spo_bou_ons] = ...
 noi_sig, ext_eth_dat(sta_tim, end_tim, fra_tim, fra, fra_rat_get, ...
                mas_hea, mas_gil, mas_mou, mas_eye, mas_noi, ...
                tai_bea, pri_com_hea, pri_com_ope, pri_com_mou, pri_com_eye, ...
                hea_per_pro, ope_per_pro, mou_per_pro, ...
                ope_per_wid_min, ope_per_dec, ope_per_wid_max, ...
                ope_dir, mou_dir, min_thr, ...
                dur_spo, n_con, dur_con, isi, dur_con_end);
toc()

% =========================
%  6) PLOT: SUBPLOTS VERSION
%  =========================
dru_ons = []; sti_ons = [];
rat_tim     = tim_bin;
rat_tim_par = rat_bin_par;
tai_sig     = [];  % no continuous tail trace plotted here

plo_eth_subplots(sta_tim, end_tim, ...
    dru_ons, sti_ons, tai_bea, hea_bea, ope_bea, mou_bea, eye_bea, ...
    tai_sig, hea_sig, ope_sig, mou_sig, eye_sig, ...
    rat_tim, rat_tim_par, fra_tim, tim_rob, rat_rob, spo_bou_ons, ...
    col_tai, col_hea, col_ope, col_mou, col_eye, ...
    col_tai_sig, col_hea_sig, col_ope_sig, col_mou_sig, col_eye_sig, ...
    fon_siz, col_rob, col_dru, col_sti, col_spo_bou_ons, ...
    'Behavior summary low activity');
%% =========================
%  7) CALCULATE FREQUENCIES

sig_row = horzcat(hea_sig,ope_sig,mou_sig,eye_sig);
sig_row = sig_row(35*120:end-100*120,:);

cfg.num_lags         = 1500;
cfg.windowSize       = 1024;
cfg.maxpx            = 513;
cfg.minPeakThreshold = 1;
cfg.freqRange        = [8 15];

[acs, domf, Pxxs, Fss, isoinds] = calcsort_autocorr_freq_analysis_v3( ...
sig_row', cfg.num_lags, 120, cfg.windowSize, cfg.maxpx, cfg.minPeakThreshold, cfg.freqRange);

fprintf('%s | %s: dominant f = %.3f Hz\n', fish, cond, domf);

%% =========================
%  7) SAVE ONLY WHAT YOU NEED
%  =========================
outFile = fullfile(baseFolder, sprintf('ethogram_%s.mat', datestr(now,'yyyymmdd_HHMMSS')));
save(outFile, ...
     'hea_bea','hea_sig','ope_bea','ope_sig','mou_bea','mou_sig','eye_bea','eye_sig', ...
     'tim_bin','rat_bin_par','rat_raw_bin_par','tim_rob','rat_rob','spo_bou_ons', ...
     'mas_hea','mas_gil','mas_mou','mas_eye','mas_noi', ...
     'fra_rat_get','fra_tim','sta_tim','end_tim','rat_tim','rat_tim_par', ...
     'tai_bea','vir_spe_uni_frames','tai_ang_uni_frames','bou_ons_bou','bou_off_bou',...
     'tai_ang_uni','tai_ang','fra_tim_uni');
fprintf('Saved ethogram outputs to:\n  %s\n', outFile);

% Example:
