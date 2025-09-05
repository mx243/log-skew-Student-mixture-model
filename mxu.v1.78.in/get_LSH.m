function [lifetime, survival, hazardrt] = get_LSH(log_lifetime, plotx, usealpha)
% Usage: [lifetime, survival, hazardrt] = get_LSH(log_lifetime, plotx, usealpha)
% This function takes in a 2D log_lifetime matrix, of which each row is a
% set of log of lifetime pdf at plotx, and return lifetime, survival,
% hazardrt.
%
% In the case usealpha ~= 0, the input lifetime will actually be the
% survival data, hence survival and hazardrt are meaningless and are
% set to [].

% Change Log:
%
%     1.1          10:oct:24    mx243      First version.
%     1.56         17:oct:24    mx243      Added usealpha.
%     1.71         25:jan:25    rfs34      Comments only changed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.get_LSH.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

lifetime = exp(log_lifetime);
if ~usealpha
    survival = 1 - cumtrapz(plotx, lifetime, 2); % Use numerical integration which is much faster than doing Monte-Carlo integration over alpha.
    hazardrt = - gradient(survival, plotx, 1) ./ survival; % Numerical differentiation.
else
    survival = [];
    hazardrt = [];
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
