function test_sample_proGamma(a, b, D, typ, N)
% Usage: test_sample_proGamma(a, b, D, typ, N)
% This function draws N samples using 'sample_proGamma.m', plot a histogram
% and compare it with the graph of f(m), where m ~ proGamma(a, b, D, typ) 
% and f is proportional to P(m).

% e.g. >> test_sample_proGamma(1, 20, 3, 2, 1000);

% Change Log:
%
%     1.1          05:jul:24    mx243      First version.
%     1.2          08:jul:24    mx243      Using hist() now.
%     1.14         26:aug:24    mx243      Added D.
%     1.16         03:sep:24    rfs34      Rejigged to normalise properly.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.test_sample_proGamma.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

% Obtain samples of m.
samples = nan(1, N);
for cnt = 1 : N
    samples(cnt) = sample_proGamma(a, b, D, typ);
end
max_m = max(samples);

Nbins = max(2, ceil(N / 10));

deltabin = (max_m + 1) / (Nbins - 1);
bins = [0 : Nbins - 1] * deltabin;

clf;
hold on;
hist(samples(:), bins);

% Plot f against m.
deltam = max_m / 1000;
if(typ == 1)
    x = 1 : deltam : 1.5 * max_m + 1;
else
    x = 0 : deltam : 1.5 * max_m + 1;
end
ly = log_density_proGamma(x, a, b, D, typ); 
ly = ly - max(ly);
y = exp(ly);
y = y / (deltam * sum(y));

plot(x, y * deltabin * N, 'Color', 'g');

hold off;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
