function y = log_density_truncated_Normal(x, mu, S, d)
% Usage: y = log_density_truncated_Normal(x, mu, S, d)
% This function returns log(P(x)) where x ~ N(mu, S)[x > d] 1D truncated Normal.

% Change Log:
%
%     1.1          1:sept:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_density_truncated_Normal.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

norm_const = (1 / 2) * (1 - erf(sqrt(S / 2) * (d - mu)));
y = (1 / 2) * (log(S) - log(2 * pi) - (x - mu) .^ 2 * S) - log(norm_const);
y(x < d) = -Inf;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
