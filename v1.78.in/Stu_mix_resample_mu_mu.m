function mu_mu = Stu_mix_resample_mu_mu(Priors, mu, C, S_mu)
% Usage: mu_mu = Stu_mix_resample_mu_mu(Priors, mu, C, S_mu)
% This function resamples mu_mu in the Stu mixture model.

% Change Log:
%
%     1.1          1:sept:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_mu_mu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

T = Priors.S_mu_mu + C * S_mu;  
tau = T \ (Priors.S_mu_mu * Priors.mu_mu_mu + S_mu * sum(mu, 2));
mu_mu = sample_Normal(tau, T);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
