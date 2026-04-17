function [domf, bp] = domf_bp(x, fs, band)
% DOMF_BP
% Returns dominant frequency and bandpower within a frequency band.

x = double(x(:)');

if all(~isfinite(x)) || numel(x) < 8
    domf = nan;
    bp   = nan;
    return
end

x = fillmissing(x, 'constant', 0);

[F, P] = pwelch(x, [], [], [], fs);

m = F >= band(1) & F <= band(2);
if ~any(m)
    domf = nan;
    bp   = nan;
    return
end

Fm = F(m);
Pm = P(m);

[~, ix] = max(Pm);
domf = Fm(ix);
bp   = trapz(Fm, Pm);
end
function [domf, bp] = domf_bp(x, fs, band)
% DOMF_BP
% Returns dominant frequency and bandpower within a frequency band.

x = double(x(:)');

if all(~isfinite(x)) || numel(x) < 8
    domf = nan;
    bp   = nan;
    return
end

x = fillmissing(x, 'constant', 0);

[F, P] = pwelch(x, [], [], [], fs);

m = F >= band(1) & F <= band(2);
if ~any(m)
    domf = nan;
    bp   = nan;
    return
end

Fm = F(m);
Pm = P(m);

[~, ix] = max(Pm);
domf = Fm(ix);
bp   = trapz(Fm, Pm);
end
