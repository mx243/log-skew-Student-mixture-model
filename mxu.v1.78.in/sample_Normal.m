function v = sample_Normal(mu, S, K)
% Usage: v = sample_Normal(mu, S, K)
% This function returns K samples from v ~ N(mu, S) as a D * K matrix.

% Change Log:
%
%     1.1          05:jul:24    mx243      First version.
%     1.2          08:jul:24    mx243      Dimension changed to D.
%     1.7          16:jul:24    mx243      Now can draw multiple samples at
%                                          once.
%     1.8          17:jul:24    mx243      Avoided repmat.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_Normal.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if(nargin == 2)
    K = 1;
elseif(nargin < 3)
    error('sample_Normal called with wrong number of arguments');
end

sz = size(mu);
D = sz(1, 1);
u = randn([D, K]);
mu = mu(:, ones(1, K));
v = forcechol(S) \ u + mu;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
