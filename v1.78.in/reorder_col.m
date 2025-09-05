function [c, n_c, ks_at_c, nonempty_c_num, p, m, mu, S, nu] = reorder_col(C, K, c, p, m, mu, S, nu)
% Usage: [c, n_c, ks_at_c, nonempty_c_num, p, m, mu, S, nu] = reorder_col(C, K, c, p, m, mu, S, nu)
% This function moves all the non-empty clusters to the front of 1, ..., C
% and reorder p, m, mu, S, nu accordingly. It also returns n_c recording 
% |{k: c_k = c}| for each c, ks_at_c a K * C matrix st. 
% ks_at_c(1 : n_c(c), c) is a col vector specifying all the k's st. c(k) = c. 

% Change Log:
%
%     1.1          27:aug:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.reorder_col.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

[bins, centers] = hist(c, 1 : C);

n_c = bins(bins > 0);
nonempty_c = centers(bins > 0);
nonempty_c_num = size(n_c, 2); % Number of non-empty clusters.

map_col = NaN(1, C);
for pos = 1 : nonempty_c_num
        map_col(nonempty_c(pos)) = pos; % Used to change c.
end

ks_at_c = NaN(K, nonempty_c_num);
cnt = zeros(1, nonempty_c_num);
for k = 1 : K
    c(k) = map_col(c(k)); % Change c.
    cnt(c(k)) = cnt(c(k)) + 1;
    ks_at_c(cnt(c(k)), c(k)) = k;
end

% Reorder p, m, mu, S, nu
nonempty_p = p(bins > 0);
empty_p = p(bins == 0);
p(1 : nonempty_c_num) = nonempty_p;
p(nonempty_c_num + 1 : C) = empty_p;

nonempty_m = m(bins > 0);
empty_m = m(bins == 0);
m(1 : nonempty_c_num) = nonempty_m;
m(nonempty_c_num + 1 : C) = empty_m;

nonempty_mu = mu(:, bins > 0);
empty_mu = mu(:, bins == 0);
mu(:, 1 : nonempty_c_num) = nonempty_mu;
mu(:, nonempty_c_num + 1 : C) = empty_mu;

nonempty_S = S(:, :, bins > 0);
empty_S = S(:, :, bins == 0);
S(:, :, 1 : nonempty_c_num) = nonempty_S;
S(:, :, nonempty_c_num + 1 : C) = empty_S;

nonempty_nu = nu(:, bins > 0);
empty_nu = nu(:, bins == 0);
nu(:, 1 : nonempty_c_num) = nonempty_nu;
nu(:, nonempty_c_num + 1 : C) = empty_nu;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
