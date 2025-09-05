function x = sample_Gamma(m, r, sz)
% Usage: x = sample_Gamma(m, r, sz)
% This function returns a matrix, whose dimension is given by sz,
% of samples from Gamma(m, r) parameterised as in the notes.

% Change Log:
%
%     1.1          08:jul:24    mx243      First version.
%     1.2          09:jul:24    rfs34      Fixed 1/r to 1./r .

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_Gamma.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if(nargin == 2)
    x = gamrnd(m, 1 ./ r);
elseif(nargin == 3)
    x = gamrnd(m, 1 ./ r, sz);
else
    error('sample_Gamma called with wrong number of arguments');
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
