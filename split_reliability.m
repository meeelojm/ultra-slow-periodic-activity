function [r_diag, r_null, perFishVals] = split_reliability(X1, X2, nsize)
% X1, X2: [neurons x time] for S1 and S8 (concatenated across fish)
% nsize : [nFish x 1] number of neurons per fish, in the same order used to build X1/X2
% Returns:
%   r_diag      : per-neuron Pearson r between matching rows (S1 vs S8)
%   r_null      : per-neuron r with a within-fish shuffled pairing (null)
%   perFishVals : cell array of per-fish r vectors

    % Align time length
    T = min(size(X1,2), size(X2,2));
    X1 = X1(:,1:T); X2 = X2(:,1:T);

    % Remove bad/flat neurons (zero variance in either split)
    good = std(X1,0,2)>0 & std(X2,0,2)>0;
    X1 = X1(good,:); X2 = X2(good,:);

    % Z-score across time (per neuron)
    X1z = zscore(X1, 0, 2);
    X2z = zscore(X2, 0, 2);

    % Row-wise Pearson correlation between matching rows (since zscored with N-1)
    r_diag = sum(X1z.*X2z, 2) ./ (T-1);

    % Build a within-fish shuffled null: permute neuron labels of X2 inside each fish
    idxStart = [1; cumsum(nsize(:))+1]; idxStart(end) = [];  % starts per fish in full list (before mask)
    idxEnd   = cumsum(nsize(:));                             % ends per fish in full list
    keepIdx  = find(good);                                   % mapping from filtered to original indices

    r_null = nan(size(r_diag));
    perFishVals = {};

    ptr = 1;
    for f = 1:numel(nsize)
        fr = idxStart(f):idxEnd(f);           % original indices for this fish
        fr_keep = intersect(fr, keepIdx);     % those that survived "good" filter
        if isempty(fr_keep), perFishVals{f} = []; %#ok<AGROW>
            continue;
        end

        % Map into the filtered matrices' row indices
        % (positions of fr_keep within keepIdx)
        [~, pos] = ismember(fr_keep, keepIdx);
        pos(pos==0) = [];  % safety

        % Shuffle pairing inside the fish
        shuf = pos(randperm(numel(pos)));

        % Null corr for this fish
        r_null_local = sum(X1z(pos,:).*X2z(shuf,:), 2) ./ (T-1);
        r_null(pos) = r_null_local;

        % Store per-fish real r values
        perFishVals{f} = r_diag(pos); %#ok<AGROW>

        ptr = ptr + numel(pos);
    end
end
