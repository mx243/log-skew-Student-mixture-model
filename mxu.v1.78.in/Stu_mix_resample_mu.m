function mu = Stu_mix_resample_mu(x, alpha, n_c, ks_at_c, nonempty_c_num, S, nu, S_mu, mu_mu, D, C, t_SA)
% Usage: mu = Stu_mix_resample_mu(x, alpha, n_c, ks_at_c, nonempty_c_num, S, nu, S_mu, mu_mu, D, C, t_SA)
% This function resamples phi.mu in the Stu mixture model.
% nonempty_c_num is the number of non-empty clusters.
% n_c is of size [1, nonempty_c_num] and contains the number of patients in each cluster.
% ks_at_c has size [K, nonempty_c_num], where K is the number of patients;
%  it contains the list of patients in each cluster.

% Change Log:
%
%     1.1          31:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Parameter descriptions added; ksi changed to xi to conform
%                                          with usual English and LaTeX spelling of \xi; half as 
%                                          many sqrts called.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_mu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

mu = NaN(D, C);

for it_c = 1 : C
    if it_c <= nonempty_c_num
        alpha_k = alpha(ks_at_c(1 : n_c(it_c), it_c)');
        sqrt_alpha_k = sqrt(alpha_k);
        nu_c = nu(:, it_c);
        T = S_mu + t_SA * sum(alpha_k, 2) * S(:, :, it_c);
        xi = sum(alpha_k(ones(D, 1), :) .* x(:, ks_at_c(1 : n_c(it_c), it_c)') - ...
                  sqrt_alpha_k(ones(D, 1), :) .* nu_c(:, ones(1, n_c(it_c))), 2);
        xi = T \ (S_mu * mu_mu + t_SA * (S(:, :, it_c) * xi));
        mu(:, it_c) = sample_Normal(xi, T);
    else
        mu(:, it_c) = sample_Normal(mu_mu, S_mu);
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
