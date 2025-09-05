function R = forcechol(A)
% Usage: R = forcechol(A)
% This function computes the Cholesky decomposition of a positive-definite 
% symmetric matrix A, using approximation A + eps * I to avoid
% ill-conditioned input. Also, we forced slightly asymmetric A to be
% symmetric by setting A = (A + A') / 2.

% Change Log:
%
%     1.1          19:jun:24    mx243    First version.
%     1.2          19:jun:24    rfs      Layout only adjusted.
%     1.3          24:jun:24    mx243    Treatment for zero/empty input.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.forcechol.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

% Treatment for zero/empty input
if(isempty(A) | A == zeros(size(A)))
    % R = zeros(size(A));
    error('forcechol called with zero/empty matrix');
    return;
end

A = (A + A') / 2;

eps = min(abs(diag(A)));

% Forced Cholesky
[R, flag] = chol(A);
while(flag)
    A = A + eps * eye(size(A));
    eps = eps * 2;
    [R, flag] = chol(A);
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
