function S = Stu_mix_resample_S(Priors, x, alpha, n_c, ks_at_c, nonempty_c_num, ...
                                mu, nu, D, C, m_S, R_S, t_SA)
% Usage: S = Stu_mix_resample_S(Priors, x, alpha, n_c, ks_at_c, nonempty_c_num, ...
%                               mu, nu, D, C, m_S, R_S, t_SA)
% This function resamples phi.S in the Stu mixture model.
% It assumes that non-empty clusters are ordered first.
% nonempty_c_num is the number of non-empty clusters.
% n_c is of size [1, nonempty_c_num] and contains the number of patients in each cluster.
% ks_at_c has size [K, nonempty_c_num], where K is the number of patients;
%  it contains the list of patients in each cluster.

% Change Log:
%
%     1.1          31:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Comments only changed.
%     1.18         05:sep:24    mx243      Now only resamples the non-empty
%                                          clusters, but does resize S according to the new C.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_S.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

S = NaN(D, D, C);

for it_c = 1 : C
    if it_c <= nonempty_c_num
        mu_c = mu(:, it_c);
        nu_c = nu(:, it_c);
        sqrt_alpha = sqrt(alpha(ks_at_c(1 : n_c(it_c), it_c)'));
        tmp_mat = x(:, ks_at_c(1 : n_c(it_c), it_c)') - mu_c(:, ones(1, n_c(it_c))) ...
                                                      - (nu_c(:, ones(1, n_c(it_c))) ./ sqrt_alpha(ones(D, 1), :));
        half_talpha = (t_SA / 2) * alpha(ks_at_c(1 : n_c(it_c), it_c)');
        T = nu(:, it_c) * (nu(:, it_c)' / (2 * Priors.kappa_nu)) + ...
            (tmp_mat .* half_talpha(ones(D, 1), :)) * tmp_mat';
        S(:, :, it_c) = sample_Wis(m_S + (t_SA * n_c(it_c) + 1) / 2, (m_S - 1) * R_S + T, D);
    elseif 0
        T = nu(:, it_c) * nu(:, it_c)' / (2 * Priors.kappa_nu);
        S(:, :, it_c) = sample_Wis(m_S + 1 / 2, (m_S - 1) * R_S + T, D);
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
