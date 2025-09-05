function C = sample_discrete_Exp(kappa, sz, max_C)
% Usage: C = sample_discrete_Exp(kappa, sz, max_C)
% This function returns a 1 * sz vector of samples from the discrete distr
% Exp(kappa) starting at 1 truncated above at max_C.

% Change Log:
%
%     1.1          26:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Comments only changed.
%     1.30         16:sep:24    rfs34      max_C added.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_discrete_Exp.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if (nargin < 2 || isempty(sz))
    sz = 1;
end
if (nargin < 3 || isempty(max_C)),
    max_C = Inf;
end

if nargin < 1,   
    error('sample_discrete_Exp called with wrong number of arguments');
end

tmp = rand(1, sz);
C = ceil(log(1 - tmp * (1 - kappa .^ max_C)) / log(kappa));

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
