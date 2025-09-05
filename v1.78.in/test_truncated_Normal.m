function test_truncated_Normal(mu, S, d, N)
% Usage: test_truncated_Normal(mu, S, d, N)

% Change Log:
%
%     1.1          27:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Rejigged to normalise properly.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.test_truncated_Normal.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

v = NaN(1, N);
for j = 1 : N
    v(j) = sample_truncated_Normal(mu, S, d);
end

lb = min(v);
ub = max(v);

Nbins = max(2, ceil(N / 10));

deltabin = (ub - lb) / Nbins;

bins = [lb - 0.2 * (ub - lb) : deltabin : ub + 0.2 * (ub - lb)];

clf;
hold on;

hist(v, bins);

deltax = (ub - lb) / 1000;
x = [lb - 0.2 * (ub - lb) : deltax : ub + 0.2 * (ub - lb)];

f = @(x) sqrt(S / (2 * pi)) * exp(- (1 / 2) * S * (x - mu) .^ 2) / normcdf(d, mu, 1 / sqrt(S), 'upper') .* (x > d);

plot(x, f(x) * N * deltabin, 'g');

hold off;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
