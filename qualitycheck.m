%%
% Define the base folder path
clear all, clc
baseFolder = 'E:\Emiliano\ccn2a\high2low\20260211_14_35_36_f026_ejm_ccn2a_gcamp6s_21dpf_ong_lowact\suite2p';

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
    resultsFile = fullfile(currentFolder, 'results.mat');
    metadataFile = fullfile(currentFolder, 'metadata_multimodal.mat');
    
    % Check if both files exist (exist returns 2 for files)
    if exist(resultsFile, 'file') == 2 && exist(metadataFile, 'file') == 2
        validFolders{end+1} = currentFolder;  %#ok<SAGROW>
    end
end

% Display the folders that contain both files
disp('Folders containing both results.mat and metadata_multimodal.mat:');
for i = 1:length(validFolders)
    disp(validFolders{i});
end
fpss = [];
for i = 1:length(validFolders)
    %i = 1:length(validFolders)
    dataFolder = validFolders{i}; 

    iniFiles = dir(fullfile(dataFolder, '*.ini'));
    if ~isempty(iniFiles)
        iniFilePath = fullfile(dataFolder, iniFiles(1).name);
        iniContent = fileread(iniFilePath);
        % Look for the pattern: volume.rate.(in.Hz) = <number>
        pattern = 'volume\.rate\.\(in\.Hz\)\s*=\s*([0-9.]+)';
        tokens = regexp(iniContent, pattern, 'tokens');
        if ~isempty(tokens)
            fps = str2double(tokens{1}{1});
            fpss = [fpss;fps];
        else
            fps = 14.64;
            fpss = [fpss;fps];
            warning('Could not find volume.rate.(in.Hz) in %s, using default fps = 14.64', dataFolder);
        end
    else
        fps = 14.64;
        fpss = [fpss;fps];
        warning('No .ini file found in %s, using default fps = 14.64', dataFolder);
    end
end
%% Load Data and Set Parameters
% Load suite2p results containing the fluorescence traces.
for i = 1:size(validFolders,2)%i=1:size(validFolders,2)
    %i=1:size(validFolders,2)
    dataFolder = validFolders{i}; 
    load(fullfile(dataFolder, "results.mat")); 
    initialVars = who;

    traces = trace;  % Rename for clarity
    fps = fpss(i);     % Current acquisition frame rate

    %get name of the experiment
    % Example folder name
    folderName = dataFolder;
    % Define a pattern that looks for an underscore, then the fish ID (f0 plus digits), then _ejm
    pattern = '_(f0\d+)_ejm';
    % Use regexp to capture the tokens
    tokens = regexp(folderName, pattern, 'tokens');
    % Check if a token was found and extract it
    if ~isempty(tokens)
        fishID = tokens{1}{1}; % e.g., 'f019'
        disp(fishID);
    else
        warning('Fish ID not found in the string.');
    end
    
    % Compute Baseline and Delta F/F (dF/F)
    % Define baseline parameters
    tau1 = 1;      
    tau2 = 19.4;   % Adjusted for a 30.9 Hz time window for baseline estimation
    
    % Compute baseline fluorescence using an adjustable baseline function.
    baseline_all = baseline_time_adjustable(tau1, tau2, traces, fps);
    
    % Calculate dF/F as percentage change
    dff = ((traces - baseline_all) ./ baseline_all) * 100;
    
    % Remove rows containing NaNs and rows that are entirely zero
    row_non_zeros = find(all(~isnan(dff), 2));
    dff_new = double(dff(row_non_zeros, :));
    dff_new = dff_new(any(dff_new, 2), :);
    % 
    % temp = zeros(size(dff_new,1),size(dff_new,2));
    % dff_new = horzcat(temp, dff_new);
    
    % Adjust Frame Indices for Current Frame Rate
    % Original acquisition parameters (from a lower frame rate)
    original_start_frame = 200;
    original_fps = 2.5;
    
    % Calculate the corresponding start time in seconds and convert to new frame index
    start_time_sec = original_start_frame / original_fps;
    start_frame = round(start_time_sec * fps);
    
    % Define end frame (4-minute window) and number of lags (2-minute window)
    end_frame = start_frame + round(240 * fps);  % 240 seconds = 4 minutes
    num_lags  = round(120 * fps);                 % 120 seconds = 2 minutes
    
    % Frequency Analysis Parameters and Computation
    windowSize = 1024;
    maxpx = 15;
    minPeakThreshold = 2;
    freqRange = [0.01, 0.03];
    
    % Analyze autocorrelation and frequency features on the selected time segment
    [acs, dominant_frequencies, Pxxs, Fss, isoinds] = calcsort_autocorr_freq_analysis_v3(double(dff_new(:, start_frame:end_frame)), num_lags, fps, windowSize, maxpx, minPeakThreshold, freqRange);
    % sgtitle(fishID)
    % exportgraphics(fig1, fullfile(dataFolder, [fishID '_heatmaps.png']))

    % Load Stimulus Metadata and Process Stimulus Timings
    load(fullfile(dataFolder, "metadata_multimodal.mat"));
    % Convert stimulus event times to frame indices (using current fps)
    stims = floor(eve_tim * fps);
    % Remove the first and last stimulus (potential false detections)
    stims = stims(2:end-1);
    stims_r = stims - 8000;
    % Calculate the stimulus duration in frames, scaling from an original 12 frames at 2.5 fps
    duration = floor((12 * fps) / 2.5); %o is 12 which as 5s at 2,5 hz must be 5 seconds
    
    len_res_plo_fra = floor(35*fps);
    % Extract Trial-Based dF/F Responses
    % Calculate neuron selectivity metrics using an external function
    % [dff_fra_cel, tri_lim_fra_tri_lim, f_fra_cel, dff_fra_trial_cell] = ext_dff_fra_cel_AO(trace', 0, 0, 7000, 24, stims, duration*2, len_res_plo_fra); %sti_ons should be stims
    % [dff_fra_cel_r, tri_lim_fra_tri_lim_r, f_fra_cel_r, dff_fra_trial_cell_r] = ext_dff_fra_cel_AO(trace', 0, 0, 7000, 24, stims_r, duration*2, len_res_plo_fra); %sti_ons should be stims    
    % % Replace NaNs with zeros (if any)
    % dff_fra_trial_cell(isnan(dff_fra_trial_cell)) = 0;
    % dff_fra_trial_cell_r(isnan(dff_fra_trial_cell_r)) = 0;
    % 
    % % Set Parameters for Response Cell Analysis
    % base_period_start = 1;          % Start frame of baseline period
    % base_period_end   = (duration*2)+1; % End frame of baseline period
    % stim_period_end   = floor(10 * fps); % End frame for stimulus period (10 seconds)
    % 
    % no_con = 3;  % Number of experimental conditions
    % 
    % % Define trial indices for each condition:
    % % For example, if trials 1-8 are 'light', 9-16 'tap', and 17-24 'light + tap'
    % light_indice = 1:8;
    % tap_indice   = 9:16;
    % lt_indice    = 17:24;
    % con_trials   = [light_indice; tap_indice; lt_indice];
    % 
    % 
    % % Calculate Responding Cells (Positive and Negative) Using dF/F Data
    % % --- Using the STD method (method 2) ---
    % method = 2;
    % alpha_factor = 2;  % STD factor (2 standard deviations)
    % [pos_resp_list, neg_resp_list, perc_resp_cells, perc_neg_resp_cells, info_list_dff_std] = ...
    %     calculate_resp_cells_AO(dff_fra_trial_cell, method, (duration*2)+1, ... %beore was duration+1
    %     base_period_start, base_period_end, stim_period_end, no_con, con_trials, alpha_factor);
    % 
    % method = 2;
    % alpha_factor = 2;  % STD factor (2 standard deviations)
    % [pos_resp_list_r, neg_resp_list_r, perc_resp_cells_r, perc_neg_resp_cells_r, info_list_dff_std_r] = ...
    %     calculate_resp_cells_AO(dff_fra_trial_cell_r, method, (duration*2)+1, ... %beore was duration+1
    %     base_period_start, base_period_end, stim_period_end, no_con, con_trials, alpha_factor);

    % % --- Using the Sign Rank method (method 1) ---
    % method = 1;
    % alpha_factor = 0.05;  % p-value threshold for the sign rank test
    % [pos_resp_list, neg_resp_list, perc_resp_cells, perc_neg_resp_cells, info_list_dff_signrank] = ...
    %     calculate_resp_cells_AO(dff_fra_trial_cell, method, duration+1, ...
    %     base_period_start, base_period_end, stim_period_end, no_con, con_trials, alpha_factor);
    
    % % Extract and Plot Modulated Responses (Mean Responses Per Condition)
    % % Identify modulated neurons (both positive and negative responses)
    % mod_light = find(info_list_dff_std(:, 1) ~= 0);
    % mod_tap   = find(info_list_dff_std(:, 2) ~= 0);
    % mod_lt    = find(info_list_dff_std(:, 3) ~= 0);
    % 
    % prop_mod_light = (size(mod_light,1)/size(dff_new,1))*100;
    % prop_mod_tap = (size(mod_tap,1)/size(dff_new,1))*100;
    % prop_mod_lt = (size(mod_lt,1)/size(dff_new,1))*100;
    % 
    % % Compute the mean response across trials for each condition
    % sorting_mean_light_dstd = squeeze(mean(dff_fra_trial_cell(:, 1:8, mod_light), 2));
    % sorting_mean_tap_dstd   = squeeze(mean(dff_fra_trial_cell(:, 9:16, mod_tap), 2));
    % sorting_mean_lt_dstd    = squeeze(mean(dff_fra_trial_cell(:, 17:24, mod_lt), 2));
    % 
    % sorting_mean_alllight_dstd = squeeze(mean(dff_fra_trial_cell(:, 1:8, :), 2));
    % sorting_mean_alltap_dstd = squeeze(mean(dff_fra_trial_cell(:, 9:16, :), 2));
    % sorting_mean_alllt_dstd = squeeze(mean(dff_fra_trial_cell(:, 17:24, :), 2));
    % 
    % sorting_mean_s1_alllight_dstd = squeeze(mean(dff_fra_trial_cell(:, 1, :), 2));
    % sorting_mean_s1_alltap_dstd = squeeze(mean(dff_fra_trial_cell(:, 9, :), 2));
    % sorting_mean_s1_alllt_dstd = squeeze(mean(dff_fra_trial_cell(:, 17, :), 2));    
    % 
    % sorting_mean_s8_alllight_dstd = squeeze(mean(dff_fra_trial_cell(:, 8, :), 2));
    % sorting_mean_s8_alltap_dstd = squeeze(mean(dff_fra_trial_cell(:, 16, :), 2));
    % sorting_mean_s8_alllt_dstd = squeeze(mean(dff_fra_trial_cell(:, 24, :), 2));    



    % % Plot the mean responses (one plot per condition)
    % figure
    % subplot(1, 3, 1)
    % plot(sorting_mean_light_dstd)
    % xline(duration, 'k', 'LineWidth', 3)
    % title('Modulated Light')
    % 
    % subplot(1, 3, 2)
    % plot(sorting_mean_tap_dstd)
    % xline(duration, 'k', 'LineWidth', 3)
    % title('Modulated Tap')
    % 
    % subplot(1, 3, 3)
    % plot(sorting_mean_lt_dstd)
    % xline(duration, 'k', 'LineWidth', 3)
    % title('Modulated Light + Tap')
    
    %% Plot Heatmaps of Positive and Negative Responses
    % Extract neurons with positive responses using the STD analysis results
    % mod_plight_dstd = find(info_list_dff_std(:, 1) == 1);
    % mod_ptap_dstd   = find(info_list_dff_std(:, 2) == 1);
    % mod_plt_dstd    = find(info_list_dff_std(:, 3) == 1);
    % 
    % prop_mod_plight = (size(mod_plight_dstd,1)/size(dff_new,1))*100;
    % prop_mod_ptap = (size(mod_ptap_dstd,1)/size(dff_new,1))*100;
    % prop_mod_plt = (size(mod_plt_dstd,1)/size(dff_new,1))*100;
    % 
    % sorting_mean_plight_dstd = squeeze(mean(dff_fra_trial_cell(:, 1:8, mod_plight_dstd), 2));
    % sorting_mean_ptap_dstd   = squeeze(mean(dff_fra_trial_cell(:, 9:16, mod_ptap_dstd), 2));
    % sorting_mean_plt_dstd    = squeeze(mean(dff_fra_trial_cell(:, 17:24, mod_plt_dstd), 2));
    % 
    % % Extract neurons with negative responses using the STD analysis results
    % mod_nlight_dstd = find(info_list_dff_std(:, 1) == -1);
    % mod_ntap_dstd   = find(info_list_dff_std(:, 2) == -1);
    % mod_nlt_dstd    = find(info_list_dff_std(:, 3) == -1);
    % 
    % prop_mod_nlight = (size(mod_nlight_dstd,1)/size(dff_new,1))*100;
    % prop_mod_ntap = (size(mod_ntap_dstd,1)/size(dff_new,1))*100;
    % prop_mod_nlt = (size(mod_nlt_dstd,1)/size(dff_new,1))*100;    
    % 
    % sorting_mean_nlight_dstd = squeeze(mean(dff_fra_trial_cell(:, 1:8, mod_nlight_dstd), 2));
    % sorting_mean_ntap_dstd   = squeeze(mean(dff_fra_trial_cell(:, 9:16, mod_ntap_dstd), 2));
    % sorting_mean_nlt_dstd    = squeeze(mean(dff_fra_trial_cell(:, 17:24, mod_nlt_dstd), 2));
    % 
    % sorting_means.p_light = sorting_mean_plight_dstd;
    % sorting_means.p_tap   = sorting_mean_ptap_dstd;
    % sorting_means.p_lt    = sorting_mean_plt_dstd;
    % sorting_means.n_light = sorting_mean_nlight_dstd;
    % sorting_means.n_tap   = sorting_mean_ntap_dstd;
    % sorting_means.n_lt    = sorting_mean_nlt_dstd;

%% plot

    % Create heatmaps for each condition (positive and negative responses)
    % fig1 = figure('Renderer','painters','Position',[765, 420, 765, 420]);  % Capture the figure handle
    % subplot(2, 3, 1)
    % imagesc(sorting_mean_plight_dstd')
    % xline(duration*2, 'w', 'LineWidth', 3)
    % title(['Positive Light ' num2str(prop_mod_plight) '%'])
    % colormap("turbo"); colorbar; caxis([-5 50])
    % 
    % subplot(2, 3, 2)
    % imagesc(sorting_mean_ptap_dstd')
    % xline(duration, 'w', 'LineWidth', 3)
    % title(['Positive Tap ' num2str(prop_mod_ptap) '%'])
    % colormap("turbo"); colorbar; caxis([-5 50])
    % 
    % subplot(2, 3, 3)
    % imagesc(sorting_mean_plt_dstd')
    % xline(duration, 'w', 'LineWidth', 3)
    % title(['Positive Light + Tap' num2str(prop_mod_plt) '%'])
    % colormap("turbo"); colorbar; caxis([-5 50])
    % 
    % subplot(2, 3, 4)
    % imagesc(sorting_mean_nlight_dstd')
    % xline(duration, 'w', 'LineWidth', 3)
    % title(['Negative Light ' num2str(prop_mod_nlight) '%'])
    % colormap("turbo"); colorbar; caxis([-5 50])
    % 
    % subplot(2, 3, 5)
    % imagesc(sorting_mean_ntap_dstd')
    % xline(duration, 'w', 'LineWidth', 3)
    % title(['Negative Tap ' num2str(prop_mod_ntap) '%'])
    % colormap("turbo"); colorbar; caxis([-5 50])
    % 
    % subplot(2, 3, 6)
    % imagesc(sorting_mean_nlt_dstd')
    % xline(duration, 'w', 'LineWidth', 3)
    % title(['Negative Light + Tap ' num2str(prop_mod_nlt) '%'])
    % colormap("turbo"); colorbar; caxis([-5 50])
    % 
    % sgtitle(fishID)

    % exportgraphics(fig1, fullfile(dataFolder, [fishID '_heatmaps_extended.png']))
    
    %% --- Subplot for Modulated Light ---
    %fig2 = figure('Renderer', 'painters', 'Position', [765, 420, 765, 420]);
    % 
    % % ========== 1) LIGHT SUBPLOT ==========
    % subplot(1, 3, 1); hold on
    % % Determine number of positive (plight) and negative (nlight) traces
    % nPlight = size(sorting_mean_plight_dstd, 2);
    % nNlight = size(sorting_mean_nlight_dstd, 2);
    % 
    % % ----- GREENS (Positive Light) -----
    % if nPlight > 0
    %     % Base hue around 0.33 (green), vary by ±0.07
    %     hues_g = 0.33 + linspace(-0.07, 0.07, nPlight)';
    %     % Clamp hue range to [0, 1]
    %     hues_g(hues_g < 0) = 0;  
    %     hues_g(hues_g > 1) = 1;
    % 
    %     % Saturation from 0.3 to 1, value from 0.6 to 1
    %     saturations_g = linspace(0.3, 1, nPlight)';
    %     values_g      = linspace(0.6, 1, nPlight)';
    %     greens = hsv2rgb([hues_g, saturations_g, values_g]);
    % 
    %     % Plot each positive response trace in a unique green hue
    %     for j = 1:nPlight
    %         plot(sorting_mean_plight_dstd(:, j), 'Color', greens(j, :), 'LineWidth', 2)
    %     end
    % else
    %     disp('No positive light neurons to plot.');
    % end
    % 
    % % ----- MAGENTAS (Negative Light) -----
    % if nNlight > 0
    %     % Base hue around 0.83 (magenta), vary by ±0.1
    %     hues_m = 0.83 + linspace(-0.1, 0.1, nNlight)';
    %     % Clamp hue range to [0, 1]
    %     hues_m(hues_m < 0) = 0;
    %     hues_m(hues_m > 1) = 1;
    % 
    %     % Saturation from 0.2 to 1, value from 0.5 to 1
    %     saturations_m = linspace(0.2, 1, nNlight)';
    %     values_m      = linspace(0.5, 1, nNlight)';
    %     magentas = hsv2rgb([hues_m, saturations_m, values_m]);
    % 
    %     % Plot each negative response trace in a unique magenta hue
    %     for j = 1:nNlight
    %         plot(sorting_mean_nlight_dstd(:, j), 'Color', magentas(j, :), 'LineWidth', 2)
    %     end
    % else
    %     disp('No negative light neurons to plot.');
    % end
    % 
    % xline(duration, 'k', 'LineWidth', 3)
    % xlabel('Frames'); ylabel('% \Delta F/F')
    % ylim([-50 100]); box off
    % title(['Modulated Light ' num2str(prop_mod_light) '%'])
    % 
    % % ========== 2) TAP SUBPLOT ==========
    % subplot(1, 3, 2); hold on
    % % Determine number of positive (ptap) and negative (ntap) traces
    % nPtap = size(sorting_mean_ptap_dstd, 2);
    % nNtap = size(sorting_mean_ntap_dstd, 2);
    % 
    % % ----- GREENS (Positive Tap) -----
    % if nPtap > 0
    %     hues_g = 0.33 + linspace(-0.1, 0.1, nPtap)';   % vary hue more
    %     hues_g(hues_g < 0) = 0;
    %     hues_g(hues_g > 1) = 1;
    %     saturations_g = linspace(0.2, 1, nPtap)';
    %     values_g      = linspace(0.5, 1, nPtap)';
    %     greens = hsv2rgb([hues_g, saturations_g, values_g]);
    % 
    %     for j = 1:nPtap
    %         plot(sorting_mean_ptap_dstd(:, j), 'Color', greens(j, :), 'LineWidth', 2)
    %     end
    % else
    %     disp('No positive tap neurons to plot.');
    % end
    % 
    % % ----- MAGENTAS (Negative Tap) -----
    % if nNtap > 0
    %     hues_m = 0.83 + linspace(-0.1, 0.1, nNtap)';
    %     hues_m(hues_m < 0) = 0;
    %     hues_m(hues_m > 1) = 1;
    %     saturations_m = linspace(0.2, 1, nNtap)';
    %     values_m       = linspace(0.5, 1, nNtap)';
    %     magentas = hsv2rgb([hues_m, saturations_m, values_m]);
    % 
    %     for j = 1:nNtap
    %         plot(sorting_mean_ntap_dstd(:, j), 'Color', magentas(j, :), 'LineWidth', 2)
    %     end
    % else
    %     disp('No negative tap neurons to plot.');
    % end
    % 
    % xline(duration, 'k', 'LineWidth', 3)
    % xlabel('Frames'); ylabel('% \Delta F/F')
    % ylim([-50 100]); box off
    % title(['Modulated Tap ' num2str(prop_mod_tap) '%'])
    % 
    % % ========== 3) LIGHT+TAP SUBPLOT ==========
    % subplot(1, 3, 3); hold on
    % % Determine number of positive (plt) and negative (nlt) traces
    % nPlt = size(sorting_mean_plt_dstd, 2);
    % nNlt = size(sorting_mean_nlt_dstd, 2);
    % 
    % % ----- GREENS (Positive Light+Tap) -----
    % if nPlt > 0
    %     hues_g = 0.33 + linspace(-0.1, 0.1, nPlt)';
    %     hues_g(hues_g < 0) = 0;
    %     hues_g(hues_g > 1) = 1;
    %     saturations_g = linspace(0.2, 1, nPlt)';
    %     values_g      = linspace(0.5, 1, nPlt)';
    %     greens = hsv2rgb([hues_g, saturations_g, values_g]);
    % 
    %     for j = 1:nPlt
    %         plot(sorting_mean_plt_dstd(:, j), 'Color', greens(j, :), 'LineWidth', 2)
    %     end
    % else
    %     disp('No positive light+tap neurons to plot.');
    % end
    % 
    % % ----- MAGENTAS (Negative Light+Tap) -----
    % if nNlt > 0
    %     hues_m = 0.83 + linspace(-0.1, 0.1, nNlt)';
    %     hues_m(hues_m < 0) = 0;
    %     hues_m(hues_m > 1) = 1;
    %     saturations_m = linspace(0.2, 1, nNlt)';
    %     values_m       = linspace(0.5, 1, nNlt)';
    %     magentas = hsv2rgb([hues_m, saturations_m, values_m]);
    % 
    %     for j = 1:nNlt
    %         plot(sorting_mean_nlt_dstd(:, j), 'Color', magentas(j, :), 'LineWidth', 2)
    %     end
    % else
    %     disp('No negative light+tap neurons to plot.');
    % end
    % 
    % xline(duration, 'k', 'LineWidth', 3)
    % xlabel('Frames'); ylabel('% \Delta F/F')
    % ylim([-50 100]); box off
    % title(['Modulated Light + Tap ' num2str(prop_mod_lt) '%'])
    % 
    % sgtitle(fishID)
    % 
    % exportgraphics(fig2, fullfile(dataFolder, [fishID '_lineplots_extended.pdf']));
    % 


  %% Example DFQ data for each row and condition (adapt these to your actual variables):

    % Run analysis for positively modulated light group
    % if isempty(mod_plight_dstd)
    %     groupResults.plight = struct('ongoing', [], 'light', [], 'tap', [], 'light_plus_tap', [], 'ongoing_after', []);
    %     disp(['No positively modulated light neurons found for ' fishID]);
    % else
    %     groupResults.plight = analyzeGroup(dff_new, mod_plight_dstd, stims, start_time_sec, fps, ...
    %         num_lags, windowSize, maxpx, minPeakThreshold, freqRange, fishID);
    % end
    % 
    % % Run analysis for positively modulated tap group
    % if isempty(mod_ptap_dstd)
    %     groupResults.ptap = struct('ongoing', [], 'light', [], 'tap', [], 'light_plus_tap', [], 'ongoing_after', []);
    %     disp(['No positively modulated tap neurons found for ' fishID]);
    % else
    %     groupResults.ptap = analyzeGroup(dff_new, mod_ptap_dstd, stims, start_time_sec, fps, ...
    %         num_lags, windowSize, maxpx, minPeakThreshold, freqRange, fishID);
    % end
    % 
    % % Run analysis for positively modulated light+tap group
    % if isempty(mod_plt_dstd)
    %     groupResults.plt = struct('ongoing', [], 'light', [], 'tap', [], 'light_plus_tap', [], 'ongoing_after', []);
    %     disp(['No positively modulated light+tap neurons found for ' fishID]);
    % else
    %     groupResults.plt = analyzeGroup(dff_new, mod_plt_dstd, stims, start_time_sec, fps, ...
    %         num_lags, windowSize, maxpx, minPeakThreshold, freqRange, fishID);
    % end
    % 
    % % Run analysis for negatively modulated light group
    % if isempty(mod_nlight_dstd)
    %     groupResults.nlight = struct('ongoing', [], 'light', [], 'tap', [], 'light_plus_tap', [], 'ongoing_after', []);
    %     disp(['No negatively modulated light neurons found for ' fishID]);
    % else
    %     groupResults.nlight = analyzeGroup(dff_new, mod_nlight_dstd, stims, start_time_sec, fps, ...
    %         num_lags, windowSize, maxpx, minPeakThreshold, freqRange, fishID);
    % end
    % 
    % % Run analysis for negatively modulated tap group
    % if isempty(mod_ntap_dstd)
    %     groupResults.ntap = struct('ongoing', [], 'light', [], 'tap', [], 'light_plus_tap', [], 'ongoing_after', []);
    %     disp(['No negatively modulated tap neurons found for ' fishID]);
    % else
    %     groupResults.ntap = analyzeGroup(dff_new, mod_ntap_dstd, stims, start_time_sec, fps, ...
    %         num_lags, windowSize, maxpx, minPeakThreshold, freqRange, fishID);
    % end
    % 
    % % Run analysis for negatively modulated light+tap group
    % if isempty(mod_nlt_dstd)
    %     groupResults.nlt = struct('ongoing', [], 'light', [], 'tap', [], 'light_plus_tap', [], 'ongoing_after', []);
    %     disp(['No negatively modulated light+tap neurons found for ' fishID]);
    % else
    %     groupResults.nlt = analyzeGroup(dff_new, mod_nlt_dstd, stims, start_time_sec, fps, ...
    %         num_lags, windowSize, maxpx, minPeakThreshold, freqRange, fishID);
    % end
    % 
    % % Define dfqCells arrays for each group, with emptiness checks
    % 
    % conditions = {'ongoing', 'light', 'tap', 'light + tap', 'ongoing after'};
    % numCond = numel(conditions);
    % 
    % % For positively modulated responses:
    % if isempty(groupResults.plight) || ~isfield(groupResults.plight, 'ongoing') || isempty(groupResults.plight.ongoing)
    %     dfqCells_plight = {[], [], [], [], []};
    % else
    %     dfqCells_plight = {groupResults.plight.ongoing.dfq, ...
    %                        groupResults.plight.light.dfq, ...
    %                        groupResults.plight.tap.dfq, ...
    %                        groupResults.plight.light_plus_tap.dfq, ...
    %                        groupResults.plight.ongoing_after.dfq};
    % end
    % 
    % if isempty(groupResults.ptap) || ~isfield(groupResults.ptap, 'ongoing') || isempty(groupResults.ptap.ongoing)
    %     dfqCells_ptap = {[], [], [], [], []};
    % else
    %     dfqCells_ptap = {groupResults.ptap.ongoing.dfq, ...
    %                      groupResults.ptap.light.dfq, ...
    %                      groupResults.ptap.tap.dfq, ...
    %                      groupResults.ptap.light_plus_tap.dfq, ...
    %                      groupResults.ptap.ongoing_after.dfq};
    % end
    % 
    % if isempty(groupResults.plt) || ~isfield(groupResults.plt, 'ongoing') || isempty(groupResults.plt.ongoing)
    %     dfqCells_plt = {[], [], [], [], []};
    % else
    %     dfqCells_plt = {groupResults.plt.ongoing.dfq, ...
    %                     groupResults.plt.light.dfq, ...
    %                     groupResults.plt.tap.dfq, ...
    %                     groupResults.plt.light_plus_tap.dfq, ...
    %                     groupResults.plt.ongoing_after.dfq};
    % end
    % 
    % % For negatively modulated responses:
    % if isempty(groupResults.nlight) || ~isfield(groupResults.nlight, 'ongoing') || isempty(groupResults.nlight.ongoing)
    %     dfqCells_nlight = {[], [], [], [], []};
    % else
    %     dfqCells_nlight = {groupResults.nlight.ongoing.dfq, ...
    %                        groupResults.nlight.light.dfq, ...
    %                        groupResults.nlight.tap.dfq, ...
    %                        groupResults.nlight.light_plus_tap.dfq, ...
    %                        groupResults.nlight.ongoing_after.dfq};
    % end
    % 
    % if isempty(groupResults.ntap) || ~isfield(groupResults.ntap, 'ongoing') || isempty(groupResults.ntap.ongoing)
    %     dfqCells_ntap = {[], [], [], [], []};
    % else
    %     dfqCells_ntap = {groupResults.ntap.ongoing.dfq, ...
    %                      groupResults.ntap.light.dfq, ...
    %                      groupResults.ntap.tap.dfq, ...
    %                      groupResults.ntap.light_plus_tap.dfq, ...
    %                      groupResults.ntap.ongoing_after.dfq};
    % end
    % 
    % if isempty(groupResults.nlt) || ~isfield(groupResults.nlt, 'ongoing') || isempty(groupResults.nlt.ongoing)
    %     dfqCells_nlt = {[], [], [], [], []};
    % else
    %     dfqCells_nlt = {groupResults.nlt.ongoing.dfq, ...
    %                     groupResults.nlt.light.dfq, ...
    %                     groupResults.nlt.tap.dfq, ...
    %                     groupResults.nlt.light_plus_tap.dfq, ...
    %                     groupResults.nlt.ongoing_after.dfq};
    % end
    % 
    % % Tiled Layout for Positive Responses (fig3)
    % fig3 = figure('Renderer','painters','Position',[100, 100, 1100, 900]);
    % t = tiledlayout(3, 2*numCond + 1, 'TileSpacing','compact', 'Padding','compact');
    % 
    % for row = 1:3
    %     for col = 1:numCond
    % 
    %         % Compute tile index so that each condition uses 2 adjacent columns
    %         tileIndex = (row - 1)*(2*numCond + 1) + 2*(col - 1) + 1;
    %         ax = nexttile(tileIndex, [1, 2]);
    % 
    %         % Select the appropriate DFQ data and row label based on row
    %         switch row
    %             case 1
    %                 data = dfqCells_plight{col};
    %                 rowLabel = 'Frequencies Light';
    %             case 2
    %                 data = dfqCells_ptap{col};
    %                 rowLabel = 'Frequencies Tap';
    %             case 3
    %                 data = dfqCells_plt{col};
    %                 rowLabel = 'Frequencies Light + Tap';
    %         end
    % 
    %         % If the data is empty, display a placeholder text; otherwise, plot normally.
    %         if isempty(data)
    %             text(ax, 0.5, 0.5, 'No data', 'HorizontalAlignment','center','FontSize',12);
    %             axis(ax, 'off');
    %         else
    %             imagesc(ax, data);
    %             colormap(ax, 'turbo')
    %             caxis(ax, [0 0.1])  % uniform color scale for all
    %         end
    % 
    %         title(ax, conditions{col}, 'Interpreter','none')
    %         xlabel(ax, 'Neurons')
    %         ylabel(ax, rowLabel)
    %     end
    % end
    % 
    % % Add a colorbar in the last column (spanning 3 rows)
    % cbAx = nexttile(2*numCond + 1, [3, 1]);  % last column, spanning 3 rows
    % axis(cbAx, 'off')
    % colormap(cbAx, 'turbo')
    % caxis([0 0.1])
    % cb = colorbar(cbAx, 'Location','eastoutside');
    % cb.Label.String = 'DFQ';
    % 
    % sgtitle([fishID, ' positively modulated responses'])
    % exportgraphics(fig3, fullfile(dataFolder, [fishID '_positivelymodulatedfreqs.png']));

    
    % Tiled Layout for Negative Responses (fig4)
    % fig4 = figure('Renderer','painters','Position',[100, 100, 1100, 900]);
    % t = tiledlayout(3, 2*numCond + 1, 'TileSpacing','compact', 'Padding','compact');
    % 
    % for row = 1:3
    %     for col = 1:numCond
    % 
    %         tileIndex = (row - 1)*(2*numCond + 1) + 2*(col - 1) + 1;
    %         ax = nexttile(tileIndex, [1, 2]);
    % 
    %         switch row
    %             case 1
    %                 data = dfqCells_nlight{col};
    %                 rowLabel = 'Neurons Light';
    %             case 2
    %                 data = dfqCells_ntap{col};
    %                 rowLabel = 'Neurons Tap';
    %             case 3
    %                 data = dfqCells_nlt{col};
    %                 rowLabel = 'Neuronss Light + Tap';
    %         end
    % 
    %         if isempty(data)
    %             text(ax, 0.5, 0.5, 'No data', 'HorizontalAlignment','center','FontSize',12);
    %             axis(ax, 'off');
    %         else
    %             imagesc(ax, data);
    %             colormap(ax, 'turbo')
    %             caxis(ax, [0 0.1])
    %         end
    % 
    %         title(ax, conditions{col}, 'Interpreter','none')
    %         xlabel(ax, 'Frequencies')
    %         ylabel(ax, rowLabel)
    %     end
    % end
    % 
    % cbAx = nexttile(2*numCond + 1, [3, 1]);
    % axis(cbAx, 'off')
    % colormap(cbAx, 'turbo')
    % caxis([0 0.1])
    % cb = colorbar(cbAx, 'Location','eastoutside');
    % cb.Label.String = 'DFQ';
    % 
    % sgtitle([fishID, ' negatively modulated responses'])
    % exportgraphics(fig4, fullfile(dataFolder, [fishID '_negativemodulatedfreqs.pdf']));

    


    %% save basic results
    close all

    finalVars = setdiff(who, initialVars)
    save(fullfile(dataFolder, "dffs_repact_respcells.mat"))
    
end
