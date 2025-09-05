function [x, mu, accepts, tries] ...
            = Stu_mix_resample_x1_with_mu1(x, mu, D, nonempty_c_num, n_c, ks_at_c, ...
                                           is_censored, logd,  mu_0, S_0, mu_mu, S_mu, t_SA, accepts, tries)
% Usage: [x, mu, accepts, tries] ...
%           = Stu_mix_resample_x1_with_mu1(x, mu, D, nonempty_c_num, n_c, ks_at_c, ...
%                                          is_censored, logd,  mu_0, S_0, mu_mu, S_mu, t_SA, accepts, tries)
% This function is an extra step in the main MCMC palindrome that moves x(1, k)
% and mu(1, c) together if x_k is censored for all k st, c(k) = c.

% Change Log:
%
%     1.1          21:oct:24    mx243      First version.
%     1.58         24:oct:24    mx243      Dealt with 1D case.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_x1_with_mu1.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

for it_c = 1 : nonempty_c_num
    if all(is_censored(ks_at_c(1 : n_c(it_c), it_c)'))

        xk1 = x(1, ks_at_c(1 : n_c(it_c), it_c)');

        lb = mu(1, it_c) + max(logd(ks_at_c(1 : n_c(it_c), it_c)') - xk1); % Lower bound of mu_c_1_prime.

        % Arguments of the Normal distr P(mu_c_1_prime | mu_mu, S_mu, mu(it_c)_1bar).
        if D > 1
            mu2 = mu_mu(1) - S_mu(1, 2 : D) * (mu(2 : D, it_c) - mu_mu(2 : D)) / S_mu(1, 1);
        else
            mu2 = mu_mu(1);
        end
        S2 = S_mu(1, 1);

        mu_c_1_prime = sample_truncated_Normal(mu2, S2, lb);
        xk1_prime = xk1 + mu_c_1_prime - mu(1, it_c);
        
        % log of Hastings ratio.
        logHR = - (1 - t_SA) * (S_0 / 2) * ((xk1_prime - xk1)' * (xk1_prime + xk1 - 2 * mu_0));
        
        u = log(rand);
        tries = tries + 1;
        if u < logHR
            mu(1, it_c) = mu_c_1_prime;
            x(1, ks_at_c(1 : n_c(it_c), it_c)') = xk1_prime;
            accepts = accepts + 1;
        end
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
