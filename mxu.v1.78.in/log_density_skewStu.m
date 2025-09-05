function v = log_density_skewStu(x, mu, S, m, r, nu, D)
% Usage: v = log_density_skewStu(x, mu, S, m, r, nu, D)
% This function returns the log of the density of 
% skewStu(mu, S, m, r, nu, D) at x, where the input x can be a D * K matrix
% and the output will be a 1 * K row vector.

% Change Log:
%
%     1.1          24:sep:24    mx243      First version.
%     1.45         03:oct:24    rfs34      Layout only changed.
%     1.50         07:oct:24    mx243      Return -Inf when any entry of x is Inf or -Inf.
%     1.51         09:oct:24    mx243      Return 0 when x = [].
%     1.52         14:oct:24    rfs34      If x is empty now returns an array of size [1, size(x, 2:end)].
%     1.53         16:oct:24    rfs34      Fixed stupid typo.
%     1.66         01:nov:24    mx243      Reduced the number of Eadd and
%                                          Esub called to increase speed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_density_skewStu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if isempty(x)
   % We return an empty array of size [1, size(x, 2 : end)].
   sz = size(x);
   v = NaN([1, sz(2 : end)]);
   return;
end

K = size(x, 2);
y = x - mu(:, ones(1, K));

tmpv1 = r + (1 / 2) * sum(y .* (S * y), 1); % tmpv1 = r + (1 / 2)y'Sy.
tmpv2 = (S * nu)' * y; % tmpv2 = y' * S * nu.

is_neg = tmpv2 < 0;

tmpv3 = (tmpv2 .^ 2) ./ (4 * tmpv1); % tmpv3 = (y' * S * nu) ^ 2 / 4(r + (1 / 2)y'Sy).

u = gammaln(m + D / 2) + loghypergeom1F1v2((m + D / 2) * ones(1, K), 1 / 2, tmpv3);

w = log(abs(tmpv2)) - log(tmpv1) / 2 + gammaln(m + (D + 1) / 2) + ...
    loghypergeom1F1v2((m + (D + 1) / 2) * ones(1, K), 3 / 2, tmpv3);

v = - (m + D / 2) * log(tmpv1) + m * log(r) - gammaln(m) - Edet(S) / 2 - (D / 2) * log(2 * pi) ...
    - (1 / 2) * (nu' * S * nu); 

% for k = 1 : K
%     if is_neg(k)
%         v(k) = v(k) - Esub(-u(k), -w(k));
%     else
%         v(k) = v(k) - Eadd(-u(k), -w(k));
%     end
% end

u_neg = u(is_neg);
u_pos = u(~is_neg);
w_neg = w(is_neg);
w_pos = w(~is_neg);

v(is_neg) = v(is_neg) - Esub(-u_neg, -w_neg);
v(~is_neg) = v(~is_neg) - Eadd(-u_pos, -w_pos);

inf_pos = any(isinf(x), 1); % Find xs that have inf entry.
v(inf_pos) = -Inf;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
