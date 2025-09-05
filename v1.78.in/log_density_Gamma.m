function y = log_density_Gamma(x, m, r)
% Usage: y = log_density_Gamma(x, m, r)
% This function returns log(P(x)) where x ~ Gamma(m, r). This
% function is used when integrating wrt x.

% Change Log:
%
%     1.1          31:jul:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_density_Gamma.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

y = m * log(r) - gammaln(m) + (m - 1) * log(x) - r * x;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
