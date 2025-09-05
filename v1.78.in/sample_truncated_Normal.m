function x = sample_truncated_Normal(mu, S, d)
% Usage: x = sample_truncated_Normal(mu, S, d)
% This function returns a sample from the 1D truncated Normal distr 
% x ~ N(mu, S)[x > d].

% Change Log:
%
%     1.1          27:aug:24    mx243      First version.
%     1.32         16:sep:24    rfs34      Now handles d = -Inf correctly.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_truncated_Normal.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if d == -Inf,
   x = randn / sqrt(S) + mu;
   return;
end

logf = @(x) - (1 / 2) * S * (x - mu) .^ 2;
dlogf = @(x) - S * (x - mu);
if d < mu
    pts = [(d + mu) / 2, mu + 1 / sqrt(S)];
else
    pts = [d + 1 / sqrt(S), d + 2 / sqrt(S)];
end
x = ars(logf, dlogf, pts, d, Inf);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
