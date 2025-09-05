function v = Edet(A)
% Usage: v = Edet(A)
% This function computes -log(det(A)) where A is a positive-definite 
% symmetric matrix. 

% Change Log:
%
%     1.1          19:jun:24    mx243    First version.
%     1.2          19:jun:24    rfs      Layout only adjusted.
%     1.3          24:jun:24    mx243    Added treatment for empty input.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Edet.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

% Treatment for empty input
if(isempty(A))
    v = 0;
    return;
end

C = forcechol(A);

v = -2 * sum(log(diag(C)));

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
