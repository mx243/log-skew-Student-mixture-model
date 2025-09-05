function nu = Stu_mix_resample_nu(Priors, x, alpha, n_c, ks_at_c, nonempty_c_num, S, mu, D, C, t_SA)
% Usage: nu = Stu_mix_resample_nu(Priors, x, alpha, n_c, ks_at_c, nonempty_c_num, S, mu, D, C, t_SA)
% This function resamples phi.nu in the Stu mixture model.

% Change Log:
%
%     1.1          31:aug:24    mx243      First version.
%     1.18         05:sep:24    mx243      Now only resamples the non-empty
%                                          clusters, but does resize nu according to the new C.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_nu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

nu = NaN(D, C);

for it_c = 1 : C
    if it_c <= nonempty_c_num
        sqrt_alpha = sqrt(alpha(ks_at_c(1 : n_c(it_c), it_c)'));
        mu_c = mu(:, it_c);
        T = (1 / Priors.kappa_nu + n_c(it_c) * t_SA) * S(:, :, it_c);
        zeta = sum(sqrt_alpha(ones(D, 1), :) .* ...
                   (x(:, ks_at_c(1 : n_c(it_c), it_c)') - mu_c(:, ones(1, n_c(it_c)))), 2);
        zeta =  (1 / (1 / (Priors.kappa_nu * t_SA) + n_c(it_c))) * zeta;
        nu(:, it_c) = sample_Normal(zeta, T);
    elseif 0
        nu(:, it_c) = sample_Normal(zeros(D, 1), S(:, :, it_c) / Priors.kappa_nu);
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
