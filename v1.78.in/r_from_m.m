function r = r_from_m(m, typ, D);
% Usage: r = r_from_m(m, typ, D);
% Sets r to 1, m - 1, or m + (D - 1) / 2 if typ is 0, 1, or 2 respectively.
% Note that m must be a row vector.
% D defaults to 1 if not passed; note that D should always be 1 if m is being 
%  used to draw a 1-d Student composer, even if the Student is D-dimensional;
%  D should match the dimension of the Wishart in the case of a Wishart, 
%  but will only affect the result if the proGamma is of type 2.

% Change Log:
%
%     1.1          13:sep:24    rfs34    First version.
%     1.52         14:oct:24    rfs34    Now allows for dimension > 1.
%     1.53         16:oct:24    rfs34    Fixed error counting parameters.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.r_from_m.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 3 || isempty(D),
   D = 1;
end

if size(m, 1) ~= 1,
   error('r_from_m expects a row vector');
end

rs_from_m = [ones(size(m)); m - 1; m + (D - 1) / 2];
r = rs_from_m(typ + 1, :);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
