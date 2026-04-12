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