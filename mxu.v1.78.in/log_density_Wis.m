function y = log_density_Wis(S, m, R)
% Usage: y = log_density_Wis(S, m, R)
% This function returns log(P(S)) where S ~ Wis(m, R, D@1) 1D Wishart.
% Not sure why we need this function - it just returns the Gamma density,
%  and it doesn't appear to be used anywhere.

% Change Log:
%
%     1.1          02:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Comments only changed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_density_Wis.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

y = m * log(R) - gammaln(m) + (m - 1) * log(S) - R * S;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
