function x = Stu_mix_resample_x1(x, alpha, D, K, mu, S, nu, c, logd, is_censored, mu_0, S_0, t_SA)
% Usage: x = Stu_mix_resample_x1(x, alpha, D, K, mu, S, nu, c, logd, is_censored, mu_0, S_0, t_SA)
% This function resamples x1 in the Stu mixture model.

% Change Log:
%
%     1.1          28:aug:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_x1.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

for k = 1 : K
    if is_censored(k)
        if D > 1
            A = (x(2 : D, k) - (mu(2 : D, c(k)) + nu(2 : D, c(k)) / sqrt(alpha(k))))' * S(2 : D, 1, c(k)) / ...
                S(1, 1, c(k));
        else
            A = 0;
        end
        T = t_SA * alpha(k) * S(1, 1, c(k)) + (1 - t_SA) * S_0;
        B = (t_SA * alpha(k) * S(1, 1, c(k)) * (mu(1, c(k)) + nu(1, c(k)) / sqrt(alpha(k)) - A) + ...
             (1 - t_SA) * S_0 * mu_0) / T;
        x(1, k) = sample_truncated_Normal(B, T, logd(k));
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
