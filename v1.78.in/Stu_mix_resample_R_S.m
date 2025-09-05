function R_S = Stu_mix_resample_R_S(Priors, D, C, S, m_S)
% Usage: R_S = Stu_mix_resample_R_S(Priors, D, C, S, m_S)
% This function resamples R_S in the Stu mixture model.

% Change Log:
%
%     1.1          1:sept:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_R_S.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

R_S = sample_Wis(Priors.m_R_S + C * (m_S + (D - 1) / 2), Priors.R_R_S + (m_S - 1) * sum(S, 3), D);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
