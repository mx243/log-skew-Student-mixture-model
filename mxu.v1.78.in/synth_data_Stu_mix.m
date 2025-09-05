function [x_ob, Args_true] = synth_data_Stu_mix(Priors, K, mu_1, S_1)
% Usage: [x_ob, Args_true] = synth_data_Stu_mix(Priors, K, mu_1, S_1)
% This function samples a theta in the Stu mix model, draw K samples of x.
% We then draw a logd from N(mu_1, S_1) for each k, and if
% 1. logd < x(1, k), set x_ob(1, k) = - exp(d), and patient k is censored
% 2. logd > x(1, k), set x_ob(1, k) = exp(x(1, k)), uncensored.

% Change Log:
%
%     1.1          02:sep:24    mx243      First version.
%     1.19         06:sep:24    mx243      Censoring time now indep to
%                                          data. Added n_c, ks_at_c, 
%                                          nonempty_c_num to Args_true.
%     1.20         09:sep:24    rfs34      x1 added to Args_true.
%     1.26         13:sep:24    rfs34      Now takes type of proGamma on m into account.
%     1.30         16:sep:24    rfs34      Now takes max_C into account.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.synth_data_Stu_mix.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

D = Priors.D;

% Axes above phi.
R_S_mu = sample_Wis(Priors.m_R_S_mu, Priors.R_R_S_mu, D);
S_mu = sample_Wis(Priors.m_S_mu, R_S_mu, D);
mu_mu = sample_Normal(Priors.mu_mu_mu, Priors.S_mu_mu);
m_S = sample_proGamma(Priors.a_m_S, Priors.b_m_S, D, 1);
R_S = sample_Wis(Priors.m_R_S, Priors.R_R_S, D);

% Sample phi.
C = sample_discrete_Exp(Priors.kappa_C, 1, Priors.max_C);
p = sample_Dir(Priors.kappa_eta / C, C); % It's valid to put a scalar instead of a vector of length C as the first entry.
m = NaN(1, C);
S = NaN(D, D, C);
nu = NaN(D, C);
for it_c = 1 : C
    m(it_c) = sample_proGamma(Priors.a_m, Priors.b_m, 1, Priors.typ_m); % N_m = 1.
    S(:, :, it_c) = sample_Wis(m_S, (m_S - 1) * R_S, D);
    nu(:, it_c) = sample_Normal(zeros(D, 1), S(:, :, it_c) / Priors.kappa_nu);
end
mu = sample_Normal(mu_mu, S_mu, C);

c = sample_finite_discrete(log(p), K);
[c, n_c, ks_at_c, nonempty_c_num, p, m, mu, S, nu] = reorder_col(C, K, c, p, m, mu, S, nu);

alpha = sample_Gamma(m(c), r_from_m(m(c), Priors.typ_m));

x_ob = NaN(D, K);
for k = 1 : K
    x_ob(:, k) = sample_Normal(mu(:, c(k)) + nu(:, c(k)) / sqrt(alpha(k)), alpha(k) * S(:, :, c(k)));
end

x1 = x_ob(1, :);

% Assign some patients as censored.
for k = 1 : K
    logd = sample_Normal(mu_1, S_1);
    if(logd < x_ob(1, k))
        x_ob(1, k) = - exp(logd);
    else
        x_ob(1, k) = exp(x_ob(1, k));
    end
end

Args_true = struct('x1', x1, 'alpha', alpha, 'c', c, 'n_c', n_c, ...
                   'ks_at_c', ks_at_c, 'nonempty_c_num', nonempty_c_num, ...
                   'phi', struct('C', C, 'p', p, 'm', m, ...
                                 'mu', mu, 'S', S, 'nu', nu), ...
                   'S_mu', S_mu, 'mu_mu', mu_mu, 'm_S', m_S, 'R_S', R_S, 'R_S_mu', R_S_mu);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
