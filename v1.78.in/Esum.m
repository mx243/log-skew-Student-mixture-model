function [y] = Esum(x, n);
% Usage: [y] = Esum(x, n);
% This function calculates
%  y = -log(sum(exp(-x), n))
% without ever exponentiating -x .

% Change Log:
%
%     1.1          04:jun:20    rfs      First version.
%     1.86         04:jun:20    rfs      Version number harmonisation.

% Change Log as part of mx243 code:
%
%     1.1          19:jun:24    mx243    First version.
%     1.66         01:nov:24    mx243    Removed find() to increase speed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Esum.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

shift = min(x, [], n);
ind = isinf(shift) | isnan(shift);
shift(ind) = 0;

sz = size(x);
sz = [sz, ones(n - length(sz))];
repper = ones(size(sz));
repper(n) = sz(n);

x = x - repmat(shift, repper);

y = -lognw(sum(exp(-x), n));

y = y + shift;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:

