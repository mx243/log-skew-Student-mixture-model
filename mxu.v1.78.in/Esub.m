function [y] = Esub(x1, x2);
% Usage: [y] = Esub(x1, x2);
% This function calculates
%  y = -log(exp(-x1) - exp(-x2))
% without ever exponentiating -x1 or -x2.
% It is required that all elements of x1 are less than or equal to the
%  corresponding elements of x2.

% Change Log (as /home/rfs/matlab/Esub.m):
%
%     1.1          04:jun:20    rfs      First version.
%     1.87         04:jun:20    rfs      Version number harmonisation.

% Change Log (as /home/rfs/ramakrishnan/software/lastact/mx243/Esub.m):
% 
%     1.1          19:sep:24    rfs      Taken from above.
%     1.41         23:sep:24    mx243    Force y = 0 when x1 > x2 but 
%                                        x1 - x2 < 1e-10.
%     1.42         29:sep:24    mx243    Fixed bug in the if condition on
%                                        line 50.
%     1.45         03:oct:24    rfs34    Adjusted handling of x1 slightly bigger than x2.
%     1.66         01:nov:24    mx243    Removed find() to increase speed.
%                                        Replaced && with & which works for arrays.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Esub.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

sz1 = size(x1);
sz2 = size(x2);

if length(sz1) ~= length(sz2),
   error('Esub called with arrays of different dimensionality');
end

if any(sz1 ~= sz2),
   error('Esub called with unequal sizes');
end

if any(x1(:) > x2(:) + 1e-10),
   error('Esub called with some of first argument greater than corresponding elements of second argument');
end

y = zeros(size(x1));

ind1 = (x1 < x2);
y(ind1) = x1(ind1) - log(1 - exp(x1(ind1) - x2(ind1)));

% If exactly equal we output -log(0).
ind3 = (x1 == x2);
y(ind3) = Inf;

% If the wrong way round, but only just, we shrink the first element to almost -log(0) but not quite,
% in order to avoid returning exactly zero in a situation in which that would cause trouble.
% (There is no really good way of handling this.)
ind4 = ((x1 <= x2 + 1e-10) & (x1 > x2));
y(ind4) = x1(ind4) + log(1e16);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:

