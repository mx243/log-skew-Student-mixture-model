function c = Stu_mix_resample_c(c, x, alpha, K, p, mu, S, nu, m, r, C, t_SA)
% Usage: c = Stu_mix_resample_c(c, x, alpha, K, p, mu, S, nu, m, r, C, t_SA)
% This function resamples c in the Stu mixture model.

% Change Log:
%
%     1.1          30:aug:24    mx243      First version.
%     1.19         09:sep:24    mx243      Redundant recalculation of devs removed.
%     1.28         13:sep:24    rfs34      Now takes account of typ_m via r.
%     1.29         16:sep:24    rfs34      Now anneals c also.
%     1.66         01:nov:24    mx243      Minor change to increase speed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_c.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

for k = 1 : K
    log_P_of_c_k = NaN(1, C);
    for it_c = 1 : C
        tmp_vec = x(:, k) - (mu(:, it_c) + nu(:, it_c) / sqrt(alpha(k)));
        S_c = S(:, :, it_c);
        log_P_of_c_k(it_c) = (t_SA / 2) * (- Edet(S_c)) - ...
                             (t_SA * alpha(k) / 2) * (tmp_vec' * S_c * tmp_vec) + ...
                             t_SA * log(p(it_c)) + m(it_c) * log(r(it_c)) - gammaln(m(it_c)) + ...
                             m(it_c) * log(alpha(k)) - r(it_c) * alpha(k);
    end
    c(k) = sample_finite_discrete(log_P_of_c_k);
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
