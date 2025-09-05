function [x, alpha] = sample_Student(mu, S, m, r, K)
% Usage: [x, alpha] = sample_Student(mu, S, m, r, K)
% This function returns K samples from x ~ Stu(mu, S, m, r) as a D * K matrix.
% It also returns the Student composer alpha in the sampling process as a
% 1 * K row vector.

% Change Log:
%
%     1.1          05:jul:24    mx243      First version.
%     1.2          08:jul:24    mx243      Corrected scale argument in gamrnd, also return alpha.
%     1.6          16:jul:24    mx243      Now can draw multiple samples at once.
%     1.8          17:jul:24    rfs34      Avoided repmat.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_Student.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if(nargin == 4)
    K = 1;
elseif(nargin < 5)
    error('sample_Student called with wrong number of arguments');
end

sz = size(S);
D = sz(1);
v = sample_Normal(zeros([D, 1]), S, K);
mu = mu(:, ones(1, K));
alpha = sample_Gamma(m, r, [1, K]);
% alpha_mat = sqrt(repmat(alpha, D, 1));
alpha_mat = sqrt(alpha(ones(D, 1), :));
x = v ./ alpha_mat + mu;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
