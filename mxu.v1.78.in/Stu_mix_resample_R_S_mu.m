function R_S_mu = Stu_mix_resample_R_S_mu(Priors, D, S_mu)
% Usage: R_S_mu = Stu_mix_resample_R_S_mu(Priors, D, S_mu)
% This function resamples R_S_mu in the Stu mixture model.

% Change Log:
%
%     1.1          1:sept:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_R_S_mu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

R_S_mu = sample_Wis(Priors.m_R_S_mu + Priors.m_S_mu + (D - 1) / 2, Priors.R_R_S_mu + S_mu, D);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
