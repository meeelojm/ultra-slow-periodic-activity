function [tai_ang, tai_ang_spe, tai_ang_uni, fra_tim_uni, vir_spe_uni, fra_tim, ...
    tai_ang_fre, ind_bor_poi, vir_spe, fra_ind, fra_rat_get, tim_tri_uni_fra, ...
    rec_tai_ang_uni_fra_sta_con, n_tri_con, n_con, tai_bea, ang_ave_con_win, ...
    nta_fra_sta_con, ave_nta_con_win, bea_lat_tri, per_pro_con, log_tri_con, bea_lat_con, ...
    ave_bas_spe, man_sec, ini_del] = ext_tai_dat_seg(fil_pat_tai, fil_pat_exp, fil_pat_ima_fis, seg_len, ...
    mph_tai, fra_mph_mpp, tim_bas_lim, dur_tai, lat_thr, win_siz)
% EXT_TAI_DAT_SEG Extracts and processes tail data from imaging and behavior files.
%
%   This function performs the following tasks:
%     1. Loads synchronization parameters and timing information.
%     2. Loads experimental metadata, including frame indices, stimulus onsets, and conditions.
%     3. Loads tail angle and speed signals from a tail file.
%     4. Computes a sum of two tail angle signals and removes baseline fluctuations.
%     5. Computes a smoothed version of the tail angle, its derivative (speed),
%        and uniform (resampled) versions.
%     6. Processes fish speed data (e.g., applying half-wave rectification).
%     7. Extracts tail beat parameters (e.g., beat times and angles).
%     8. Computes trial-by-trial and condition-averaged responses.
%     9. Computes additional metrics such as beat latency and baseline speed.
%
%   Inputs:
%     fil_pat_tai    - File path for tail data (angles, speeds)
%     fil_pat_exp    - File path for experimental parameters (metadata, stimulus onsets, etc.)
%     fil_pat_ima_fis- File path for imaging data (for synchronization)
%     seg_len        - Segment length (number of frames, e.g., 1024)
%     mph_tai        - Minimum peak height for tail beat angle detection (e.g., 0.069)
%     fra_mph_mpp    - Fraction to compute minimum peak prominence (e.g., 0.5)
%     tim_bas_lim    - Time baseline limits (e.g., [-1 0] for 1 s before stimulus)
%     dur_tai        - Duration for the tail response window (e.g., 1 s)
%     lat_thr        - Latency threshold to consider a tail flick as a response (e.g., 0.300 s)
%     win_siz        - Window size for baseline speed calculation (e.g., 120 s)
%
%   Outputs:
%     tai_ang                - Tail angle after baseline removal and smoothing
%     tai_ang_spe            - Tail angle speed (time derivative)
%     tai_ang_uni            - Uniformly sampled tail angle (full-wave rectified)
%     fra_tim_uni            - Uniform time vector for frames (after resampling)
%     vir_spe_uni            - Uniformized fish speed (processed with half-wave rectification)
%     fra_tim                - Original frame times (after cutting)
%     tai_ang_fre            - (Not used, returned as NaN)
%     ind_bor_poi            - Border point indices from baseline removal segmentation
%     vir_spe                - Fish speed at the selected frame indices
%     fra_ind                - Logical indices for frames within the desired time window
%     fra_rat_get            - Calculated frame rate from the imaging data
%     tim_tri_uni_fra        - Uniform trial time vector for tail data
%     rec_tai_ang_uni_fra_sta_con - Averaged tail angle response per condition (across trials)
%     n_tri_con              - Number of trials per condition
%     n_con                  - Number of conditions in the experiment
%     tai_bea                - Tail beat times (extracted from the tail signal)
%     ang_ave_con_win        - Average tail angle in defined time windows per condition
%     nta_fra_sta_con        - Averaged tail amplitude (or speed) per condition from trials
%     ave_nta_con_win        - Average tail amplitude in specified windows per condition
%     bea_lat_tri            - Tail beat latency per trial
%     per_pro_con            - Percentage of proper responses per condition
%     log_tri_con            - Logical matrix indicating responsive trials per condition
%     bea_lat_con            - Condition-averaged tail beat latencies
%     ave_bas_spe            - Average baseline speed (from speed data)
%     man_sec                - Manual seconds factor (from synchronization)
%     ini_del                - Initial delay (from synchronization)
%
% Written by: [Your Name], [Date]

%% --- STEP 1: Load Synchronization Parameters and Frame Timing ---
[man_sec, ini_del] = ext_syn_par(fil_pat_ima_fis, fil_pat_exp);
% Get frame indices, calculated frame rate, and frame times for the imaging data
[fra_ind, fra_rat_get, fra_tim] = ext_tem_fra_dat_cut(fil_pat_ima_fis, fil_pat_exp, man_sec, ini_del);

% Load experimental metadata (stimulus times, condition counts, etc.)
load(fil_pat_exp, 'sta_tim', 'end_tim', 'fra_rat', 'isi', 'sti_ons', 'n_con', 'n_tri_con')
n_fra = length(fra_ind);
sti_ons_tri = sti_ons;  % keep a copy of stimulus onsets for trials

%% --- STEP 2: Load Tail Data ---
load(fil_pat_tai, 'fish_angle0', 'fish_angle1', 'fish_speed')
% Use the first n_fra frames of tail angle signals
tai_ang_one = fish_angle0(1:n_fra);
tai_ang_two = fish_angle1(1:n_fra);
% Sum the two signals to get a combined tail angle signal
tai_ang_sum = tai_ang_one + tai_ang_two;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- STEP 3: Broader Frame Extraction for Synchronization ---
%[fra_ind_bro, fra_tim_bro] = ext_tem_fra_dat_cut_bro(fil_pat_ima_fis, fil_pat_exp, man_sec, ini_del);
[fra_ind_bro, fra_tim_bro] = ext_tem_fra_dat_cut_bro(fil_pat_ima_fis, fil_pat_exp, );
% Extract synchronized tail angle from the broader set of frames
tai_ang_bro = tai_ang_sum(fra_ind_bro);
% Remove the baseline dynamically using a built-in function (rem_bas_dyn)
tai_ang_bro = rem_bas_dyn(fra_rat, tai_ang_bro, seg_len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tai_ang_fre = nan;  % This output is not computed; returned as NaN

%% --- STEP 4: Baseline Removal and Smoothing for Tail Angle ---
% Get tail angle for selected frames and remove baseline
tai_ang_raw = tai_ang_sum(fra_ind);
[tai_ang, ind_bor_poi] = rem_bas_dyn(fra_rat, tai_ang_raw, seg_len);
tai_ang = smo(tai_ang);  % Smooth the baseline-corrected tail angle

% Compute the derivative (speed) of the tail angle signal
ord = 1;  % first order derivative
tai_ang_spe = dif_sig(tai_ang_raw, ord) ./ dif_sig(fra_tim, ord);  % rate of change
tai_ang_spe = abs(tai_ang_spe);  % take absolute value

%% --- STEP 5: Process Fish Speed Data ---
vir_spe_raw = fish_speed(1:n_fra);
% Select corresponding speed data from the broader frame set
vir_spe_bro = vir_spe_raw(fra_ind_bro);
% Uniformly sample and smooth the speed signal
vir_spe_uni = smo_uni(fra_tim_bro, vir_spe_bro, sta_tim, end_tim, fra_rat);
% Perform half-wave rectification on the uniform speed signal
vir_spe_uni = per_hal_wav_rec(vir_spe_uni);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thr = 0.1;  % Threshold for speed filtering
env = vir_spe_bro > thr;  % Create an envelope: 1 where speed exceeds threshold
tai_ang_bro = tai_ang_bro .* env;  % Apply envelope to tail angle

% Uniformly sample the tail angle signal (after baseline removal) from broader frames
[tai_ang_uni, fra_tim_uni] = smo_uni(fra_tim_bro, tai_ang_bro, sta_tim, end_tim, fra_rat);
% Full-wave rectify the uniform tail angle signal
tai_ang_uni = per_ful_wav_rec(tai_ang_uni);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the selected frames, extract fish speed for further analysis
vir_spe = vir_spe_raw(fra_ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
env = vir_spe > thr;  % Apply threshold to raw speed
tai_ang = tai_ang .* env;  % Multiply tail angle by the envelope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- STEP 6: Extract Tail Beat Features ---
% Extract tail beat peaks and corresponding beat angles
[tai_bea, tai_sig, tai_bea_ang] = ext_tai_bea_seg(fra_tim, mph_tai, fra_mph_mpp, tai_ang);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- STEP 7: Compute Trial-based Uniform Responses ---
% Define baseline and response durations (using inter-stimulus interval, isi)
bas_dur = isi;
res_dur = isi;
% Compute uniform trial signals (both tail angle and amplitude) for each trial
[tim_tri_uni_fra, rec_tai_ang_uni_fra_tri, n_fra_tri, nta_fra_tri] = ...
    com_rec_sig_uni_fra_tri(fra_tim_uni, tai_ang_uni, sti_ons_tri, bas_dur, res_dur, tim_bas_lim);

% Trial averaging across conditions
n_sta = 2;
rec_tai_ang_uni_fra_sta_con = nan(n_fra_tri, n_sta, n_con);
nta_fra_sta_con = nan(n_fra_tri, n_sta, n_con);
for i = 1:n_con
    % Average tail angle across trials for condition i using a built-in function
    [rec_tai_ang_uni_fra_sta_con(:, 1, i), ~, rec_tai_ang_uni_fra_sta_con(:, 2, i)] = ...
        com_sta(rec_tai_ang_uni_fra_tri(:, (i - 1)*n_tri_con + 1:i*n_tri_con), 2);
    % Average tail amplitude (or speed) across trials for condition i
    [nta_fra_sta_con(:, 1, i), ~, nta_fra_sta_con(:, 2, i)] = ...
        com_sta(nta_fra_tri(:, (i - 1)*n_tri_con + 1:i*n_tri_con), 2);
end

%% --- STEP 8: Time Averaging and Window-based Metrics ---
n_win_tot = 3;  % Number of time windows for analysis
ang_ave_con_win = nan(n_con, n_win_tot);  % Average tail angle per window per condition
ave_nta_con_win = nan(n_con, n_win_tot);   % Average tail amplitude per window per condition

% Define logical indices for different windows based on trial time vector
log_bas = tim_tri_uni_fra > -dur_tai & tim_tri_uni_fra < 0;
log_sti = tim_tri_uni_fra > 0 & tim_tri_uni_fra < dur_tai;
log_ton = tim_tri_uni_fra > dur_tai & tim_tri_uni_fra < 2*dur_tai;

for i = 1:n_con
    % Average tail angle over baseline, stimulus, and post-stimulus windows
    ang_ave_con_win(i, 1) = mean(rec_tai_ang_uni_fra_sta_con(log_bas, 1, i));
    ang_ave_con_win(i, 2) = mean(rec_tai_ang_uni_fra_sta_con(log_sti, 1, i));
    ang_ave_con_win(i, 3) = mean(rec_tai_ang_uni_fra_sta_con(log_ton, 1, i));
    %
    % Compute average tail amplitude in windows using a helper function
    ave_nta_con_win(i, :) = com_ave_amp_uni_win(tim_tri_uni_fra, ...
        reshape(nta_fra_sta_con(:, 1, i), [n_fra_tri 1]), dur_tai);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- STEP 9: Compute Beat Latency and Response Proportions ---
% Compute tail beat latency and beat angles per trial
[bea_lat_tri, ~] = ext_bea_lat_tri(sti_ons_tri, tai_bea, tai_bea_ang);
% Compute the percentage of proper responses per condition using latency threshold
[per_pro_con, log_tri_con] = com_per_pro_con(n_con, bea_lat_tri, lat_thr);

% Organize beat latencies per condition
bea_lat_con = nan(n_tri_con, n_con);
for j = 1:n_con
    bea_lat_con(:, j) = bea_lat_tri((j - 1)*n_tri_con + 1:j*n_tri_con);
end

% Remove latencies that exceed the threshold
log_nox_tri = bea_lat_tri > lat_thr;
bea_lat_tri(log_nox_tri) = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- STEP 10: Compute Average Baseline Speed ---
[~, ~, ave_bas_spe] = sel_sig_uni(fra_tim_uni, vir_spe_uni, sti_ons_tri(1), win_siz, 'sto');

end
