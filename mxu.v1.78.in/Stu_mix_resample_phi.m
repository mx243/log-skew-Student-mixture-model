function [C, p, m, mu, S, nu] = Stu_mix_resample_phi(C, p, m, mu, S, nu, ...
                                                     x, alpha, K, c, n_c, ks_at_c, nonempty_c_num, ...
                                                     Priors, S_mu, mu_mu, m_S, R_S, D, t_SA, Enables)
% Usage: [C, p, m, mu, S, nu] = Stu_mix_resample_phi(C, p, m, mu, S, nu, ...
%                                                    x, alpha, K, c, n_c, ks_at_c, nonempty_c_num, ...
%                                                    Priors, S_mu, mu_mu, m_S, R_S, D, t_SA, Enables)
% This function resamples phi in the Stu mixture model. It resamples the C
% axis using Metroplis-Hastings, with proposed distr P' ( C' | C ) = 
% \Chi(C') * 1_{floor(C / 2) <= C' <= 2 * C + 1}, where \Chi is the
% infinite-dimensional vector proportional to the vector of P(C) from which
% we're trying to sample.
% On entry the non-empty clusters must be the first in the order.
% nonempty_c_num is the number of non-empty clusters.
% n_c is of size [1, nonempty_c_num] and contains the number of patients in each cluster.
% ks_at_c has size [K, nonempty_c_num], where K is the number of patients;
%  it contains the list of patients in each cluster.

% Change Log:
%
%     1.1          30:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Parameter descriptions added; missing t_SA passed to
%                                          resample_mu; palindrome completed for non-empty clusters;
%                                          code simplified as validity of Enables.phi.C now 
%                                          checked elsewhere.
%     1.17         03:sep:24    mx243      Removed the last resampling of p, as it's independent
%                                          of the rest of phi apart from C given C and c.
%     1.18         05:sep:24    mx243      Only use Gibbs for the S_c, nu_c with non-empty c, the rest 
%                                          are resampled directly and marginally at the end.
%     1.29         16:sep:24    rfs34      Passes t_SA to resamples p, c, largeC as c now annealed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_phi.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if Enables.phi.C % This means all other axes in phi need to be resampled.
    C = Stu_mix_resample_largeC(n_c, nonempty_c_num, C, Priors, t_SA);
end

% If Enables.phi.C is set, then it is checked in sample_Stu_mix.m that the Enables.phi
% elements below must also be set.

if Enables.phi.p
    p = Stu_mix_resample_p(Priors, C, n_c, nonempty_c_num, t_SA); 
end

% Gibbs for the rest.
if Enables.phi.m
    m = Stu_mix_resample_m(Priors, alpha, n_c, ks_at_c, nonempty_c_num, C);
end

if Enables.phi.mu
    mu = Stu_mix_resample_mu(x, alpha, n_c, ks_at_c, nonempty_c_num, S, nu, S_mu, mu_mu, D, C, t_SA);
end

if Enables.phi.S
    S = Stu_mix_resample_S(Priors, x, alpha, n_c, ks_at_c, nonempty_c_num, mu, nu, D, C, m_S, R_S, t_SA);
end

if Enables.phi.nu
    nu = Stu_mix_resample_nu(Priors, x, alpha, n_c, ks_at_c, nonempty_c_num, S, mu, D, C, t_SA);
end

if Enables.phi.S
    S = Stu_mix_resample_S(Priors, x, alpha, n_c, ks_at_c, nonempty_c_num, mu, nu, D, C, m_S, R_S, t_SA);
end

if Enables.phi.mu
    mu = Stu_mix_resample_mu(x, alpha, n_c, ks_at_c, nonempty_c_num, S, nu, S_mu, mu_mu, D, C, t_SA);
end

if Enables.phi.m
    m = Stu_mix_resample_m(Priors, alpha, n_c, ks_at_c, nonempty_c_num, C);
end

[S, nu] = Stu_mix_resample_S_and_nu_empty(Priors, S, nu, m_S, R_S, nonempty_c_num, C, D, Enables);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
