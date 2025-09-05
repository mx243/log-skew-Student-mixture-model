function y = log_m_density(m, a, b)
% Usage: y = log_m_density(m, a, b)
% Returns the logarithm of the prior probability density of m,
%  given by
%  P(m) proportional to m ^ (b * m) * exp(-(a + b) * m) / Gamma(m) ^ b.

% Change Log:
%
%     1.1          28:sep:19    jc2062   As first received from JC2062.
%     1.2          10:oct:19    jc2062   As received from jc2062.
%     1.43         23:oct:19    rfs34    Comments only changed.

% Change Log (as part of mx243 code):
%
%     1.1          08:jul:24    mx243    Taken from above.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_m_density.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

y = (m .* b .* log(m)) - (b .* gammaln(m)) - (m .* (a + b));

return;

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
