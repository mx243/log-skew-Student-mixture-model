function val = F_of_theta_and_x(x, alpha, c, p, nonempty_c_num, mu, S, nu, is_censored, logd, mu_0, S_0, D, K)
% Usage: val = F_of_theta_and_x(x, alpha, c, p, nonempty_c_num, mu, S, nu, is_censored, logd, mu_0, S_0, D, K)
% This function returns the value of the expression F(\theta, x) of which 
% the expectation is the integrand in the thermodynamic integration, 
% e.g. when K = 3, patient 1 and 2 are censored while 3 is not, 
% F(\theta, x) = log( P(x(:, 1) | \theta) / P_1(x(1, 1)) *
%                     P(x(:, 2) | \theta) / P_2(x(1, 2)) *
%                     P(x(:, 3) | \theta) ),
% where P_k is the truncated Normal proportional to P_0(x)[x > logd(k)], 
% P_0 ~ N(mu_0, S_0).
% This version returns val as two elements, the first from annealing x, and the second from annealing c.

% Change Log:
%
%     1.1          01:sep:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Comments only changed.
%     1.17         03:sep:24    mx243      Added argument nonempty_c_num.
%                                          Calculate log(det(S_c_k)) first
%                                          to reduce the number of Edet called.
%     1.33         17:sep:24    rfs34      Added terms necessary when annealing c and not taking
%                                          expectation of the additional factors at the end.
%     1.33.1.2     17:sep:24    rfs34      Added p to the parameter list.
%     1.34         17:sep:24    mx243      C = length(p) on line 33.
%     1.41         23:sep:24    rfs34      Splits the returned value into the part from x1 and the part from c.
 
mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.F_of_theta_and_x.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

C = length(p);

val = NaN(2, 1);

val(1) = D * sum(log(alpha), 2) - D * K * log(2 * pi);
log_det_S_c = NaN(1, nonempty_c_num);

for it_c = 1 : nonempty_c_num % Need only to consider the non-empty clusters.
    log_det_S_c(it_c) = - Edet(S(:, :, it_c));
end

for k = 1 : K
    mu_tmp = mu(:, c(k)) + nu(:, c(k)) / sqrt(alpha(k));
    val(1) = val(1) + log_det_S_c(c(k)) ... % 1 <= c(k) <= nonempty_c_num for all k.
              - (x(:, k) - mu_tmp)' * alpha(k) * (S(:, :, c(k)) * (x(:, k) - mu_tmp));
end
val(1) = val(1) / 2;

for k = 1 : K
    if is_censored(k)
        val(1) = val(1) - log_density_truncated_Normal(x(1, k), mu_0, S_0, logd(k));
    end
end

% Terms from annealing c.
val(2) = sum(log(p(c))) + K * log(C);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
