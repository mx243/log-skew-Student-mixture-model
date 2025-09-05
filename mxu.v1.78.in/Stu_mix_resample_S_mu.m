function S_mu = Stu_mix_resample_S_mu(Priors, mu, D, C, mu_mu, R_S_mu)
% Usage: S_mu = Stu_mix_resample_S_mu(Priors, mu, D, C, mu_mu, R_S_mu)
% This function resamples S_mu in the Stu mixture model.

% Change Log:
%
%     1.1          1:sept:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_S_mu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

U = (mu - mu_mu(:, ones(1, C))) * (mu - mu_mu(:, ones(1, C)))' / 2;
S_mu = sample_Wis(Priors.m_S_mu + C / 2, R_S_mu + U, D);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
