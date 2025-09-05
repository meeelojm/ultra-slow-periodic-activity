function M = orange_map(n)
% Light peach -> deep orange
if nargin<1, n = 256; end
t = linspace(0,1,n)';
R = ones(n,1);
G = 0.95 - 0.45*t;                % 0.95 → 0.50
B = 0.85 - 0.85*t;                % 0.85 → 0.00
M = [R max(0,min(1,G)) max(0,min(1,B))];
end
