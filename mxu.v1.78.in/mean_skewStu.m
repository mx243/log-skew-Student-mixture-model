function y = mean_skewStu(mu, nu, m, r)
% Usage: y = mean_skewStu(mu, nu, m, r)
% This function returns the mean of a skew-Student distribution.

% Change Log:
%
%     1.1          18:sep:24    mx243      First version.
%     1.40         20:sep:24    rfs34      Error checking made strict.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.mean_skewStu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if m <= 1 / 2
    error('mean_skewStu called with m <= 1 / 2.');
end

y = mu + (gamma(m - (1 / 2)) / gamma(m)) * sqrt(r) * nu;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
