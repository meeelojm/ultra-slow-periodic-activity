%% Simple ethogram loader + plotter
clear; clc;

baseRoot = '\\forskning.it.ntnu.no\ntnu\mh\kin\yaksi\temp\Emilian temp\ccn2a\high2low';

fishList = {'f010','f011','f012'};
condList = {'highact','lowact'};
% so fat 3 fish out of 9 have the signals extracted, let's focus on the pipeline working and then I will extract the other signals

for iFish = 1:numel(fishList)
    fish = fishList{iFish};

    for iCond = 1:numel(condList)
        cond = condList{iCond};

        folderNow = fullfile(baseRoot, fish, cond, 'suite2p');

        if ~exist(folderNow, 'dir')
            fprintf('Missing folder: %s\n', folderNow);
            continue
        end

        files = dir(fullfile(folderNow, 'ethogram_????????_??????.mat'));;

        if isempty(files)
            fprintf('No ethogram files in: %s\n', folderNow);
            continue
        end

        for iFile = 1:numel(files)
            fileNow = fullfile(files(iFile).folder, files(iFile).name);
            fprintf('Loading: %s\n', fileNow);

            try
                load(fileNow);

                % exactly your plotting setup
                dru_ons = [];
                sti_ons = [];
                rat_tim = tim_bin;
                rat_tim_par = rat_bin_par;
                tai_sig = [];

                %figure('Color','w','Name',files(iFile).name,'NumberTitle','off');

                plo_eth_subplots(sta_tim, end_tim, ...
                    dru_ons, sti_ons, tai_bea, hea_bea, ope_bea, mou_bea, eye_bea, ...
                    tai_sig, hea_sig, ope_sig, mou_sig, eye_sig, ...
                    rat_tim, rat_tim_par, fra_tim, tim_rob, rat_rob, spo_bou_ons, ...
                    col_tai, col_hea, col_ope, col_mou, col_eye, ...
                    col_tai_sig, col_hea_sig, col_ope_sig, col_mou_sig, col_eye_sig, ...
                    fon_siz, col_rob, col_dru, col_sti, col_spo_bou_ons, ...
                    ['Behavior summary - ' fish ' - ' cond ' - ' files(iFile).name]);

            catch ME
                warning('Failed on file:\n%s\n%s', fileNow, ME.message);
            end

        end
    end
end
