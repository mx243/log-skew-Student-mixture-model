function alpha = Stu_mix_resample_alpha(alpha, x, K, mu, S, nu, m, r, c, t_SA)
% Usage: alpha = Stu_mix_resample_alpha(alpha, x, K, mu, S, nu, m, r, c, t_SA)
% This function resamples alpha in the Stu mixture model. It uses ars when
% the pdf is log-concave, otherwises it takes a Metroplis-Hastings step.

% Change Log:
%
%     1.1          28:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Parameter c added; devs now only calculated once;
%                                          unnecessary min removed.
%     1.18         04:sep:24    mx243      Changed log_alpha_k_pdf so it
%                                          works for a vector a.
%     1.19         06:sep:24    mx243      Now computes x(:, k) - mu(:, c(k)) only once.
%     1.26         13:sep:24    rfs34      Adapted to use general r as type of proGamma may not be 2.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_alpha.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

D = size(S, 1);

for k = 1 : K
    tmpv = x(:, k) - mu(:, c(k));
    C_3 = - (t_SA / 2) * nu(:, c(k))' * S(:, :, c(k)) * nu(:, c(k));
    C_2 = t_SA * D / 2 + m(c(k)) - 1;
    C_1 = (t_SA / 2) * (tmpv' * S(:, :, c(k)) * nu(:, c(k)));
    C_0 = - (t_SA / 2) * (tmpv' * S(:, :, c(k)) * tmpv) - r(c(k));

    log_alpha_k_pdf = @(a) C_2 * log(a) + 2 * C_1 * sqrt(a) + C_0 * a + C_3; % Works for vector a now.

    if (x(:, k) - mu(:, c(k)))' * S(:, :, c(k)) * nu(:, c(k)) >= 0 && t_SA * D / 2 + m(c(k)) >= 1 % ARS.

        deriv_log_alpha_k_pdf = @(a) C_2 ./ a + C_1 ./ sqrt(a) + C_0;
        
        zp = (- C_1 + sqrt(C_1 ^ 2 - 4 * C_2 * C_0)) / (2 * C_2);
        lp = 1 / (zp * 2) ^ 2;
        rp = 1 / (zp / 2) ^ 2;
        pts = [lp, rp];
        alpha(k) = ars(log_alpha_k_pdf, deriv_log_alpha_k_pdf, pts, 0, Inf);
    else % Metroplis-Hastings
        m_prop = t_SA * D / 2 + m(c(k));
        r_prop = - C_0;
        log_alpha_k_pdf_prop = @(a) (m_prop - 1) * log(a) - r_prop * a;
        alpha_prime = sample_Gamma(m_prop, r_prop);
        % The proposed distrs on a' doesn't depend on the the previous
        % state a, hence no need to normalise them when calculating the 
        % Hastings ratio (const cancels out).
        logA = log_alpha_k_pdf(alpha_prime) + log_alpha_k_pdf_prop(alpha(k)) - ...
               log_alpha_k_pdf(alpha(k)) - log_alpha_k_pdf_prop(alpha_prime);
        u = rand;
        if log(u) < logA
            alpha(k) = alpha_prime;
        end
    end
end

return;

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
