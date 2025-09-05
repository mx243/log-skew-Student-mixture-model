function [y] = Eadd(x1, x2);
% Usage: [y] = Eadd(x1, x2);
% This function calculates
%  y = -log(exp(-x1) + exp(-x2))
% without ever exponentiating either -x1 or -x2.

% Change Log:
%
%     1.1          04:jun:20    rfs      First version.
%     1.86         04:jun:20    rfs      Version number harmonisation.

% Change Log as part of mx243 code:
%
%     1.1          19:jun:24    mx243    First version.
%     1.66         01:nov:24    mx243    Removed find() to increase speed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Eadd.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

sz1 = size(x1);
sz2 = size(x2);

if length(sz1) ~= length(sz2),
   error('Eadd called with arrays of different dimensionality');
end

if any(sz1 ~= sz2),
   error('Eadd called with unequal sizes');
end

y = zeros(sz1);

ind1 = (x1 < x2);
y(ind1) = x1(ind1) - log(1 + exp(x1(ind1) - x2(ind1)));
ind2 = (x1 >= x2);
y(ind2) = x2(ind2) - log(1 + exp(x2(ind2) - x1(ind2)));
ind3 = ((x1 == x2) & isinf(x1));
y(ind3) = x1(ind3);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:

