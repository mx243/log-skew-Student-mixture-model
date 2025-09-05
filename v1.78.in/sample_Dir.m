function p = sample_Dir(eta, C)
% Usage: p = sample_Dir(eta, C)
% This function returns a 1 * C vector p which is a sample from Dir(eta, C).

% Change Log:
%
%     1.1          26:aug:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_Dir.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

p = sample_Gamma(eta, ones(1, C));
p = p / sum(p, 2);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
