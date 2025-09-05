function v = log_con_density_skewStu(x1, x1_bar, mu, S, m, r, nu, D)
% Usage: v = log_con_density_skewStu(x1, x1_bar, mu, S, m, r, nu, D)
% This function returns the log of a const multiple of the density of 
% skewStu(mu, S, m, r, nu, D) at (x1, x1_bar), where the input x1 can be a 
% vector. 
% This function can be used for numerical integration, which is then used 
% to calculate P(x1 | x1_bar).

% Change Log:
%
%     1.1          13:sep:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_con_density_skewStu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

num = length(x1);

y1 = x1 - mu(1); % y1 is a row vector because x1 is.
if D > 1
    y1_bar = x1_bar - mu(2 : D);
end

tmpv1 = r + (1 / 2) * S(1, 1) * (y1 .* y1); % tmpv1 = r + (1 / 2)y'Sy.
if D > 1
    tmpv1 = tmpv1 + (S(1, 2 : D) * y1_bar) * y1 + (1 / 2) * (y1_bar' * S(2 : D, 2 : D) * y1_bar);
end

tmpv2 = y1 * (S(1, :) * nu); % tmpv2 = y' * S * nu.
if D > 1
    tmpv2 = tmpv2 + y1_bar' * (S(2 : D, :) * nu);
end
is_neg = tmpv2 < 0;

tmpv3 = (tmpv2 .^ 2) ./ (4 * tmpv1); % tmpv3 = (y' * S * nu) ^ 2 / 4(r + (1 / 2)y'Sy).

u = gammaln(m + D / 2) + loghypergeom1F1v2((m + D / 2) * ones(1, num), 1 / 2, tmpv3);
w = log(abs(tmpv2)) - log(tmpv1) / 2 + gammaln(m + (D + 1) / 2) + ...
    loghypergeom1F1v2((m + (D + 1) / 2) * ones(1, num), 3 / 2, tmpv3);
v = - (m + D / 2) * log(tmpv1);

for j = 1 : num
    if is_neg(j)
        v(j) = v(j) - Esub(-u(j), -w(j));
    else
        v(j) = v(j) - Eadd(-u(j), -w(j));
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
