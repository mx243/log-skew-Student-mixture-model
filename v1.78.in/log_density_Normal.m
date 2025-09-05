function y = log_density_Normal(x, mu, S)
% Usage: y = log_density_Normal(x, mu, S)
% This function returns log(P(x)) where x ~ N(mu, S) 1D Gaussian. This
% function is used when integrating wrt x.

% Change Log:
%
%     1.1          31:jul:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_density_Normal.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

y = (1 / 2) * (log(S) - log(2 * pi) - ((x - mu) .^ 2) * S);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
