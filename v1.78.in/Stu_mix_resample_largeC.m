function C = Stu_mix_resample_largeC(n_c, nonempty_c_num, C, Priors, t_SA)
% Usage: C = Stu_mix_resample_largeC(n_c, nonempty_c_num, C, Priors, t_SA)
% This function resamples phi.C in the Stu mixture model. 
% It resamples C using Metroplis-Hastings, with proposed distr 
% P' ( C' | C ) = \Chi(C') * 1_{floor(C / 2) <= C' <= 2 * C + 1}, where 
% \Chi is the infinite-dimensional vector proportional to the vector of 
% P(C) from which we're trying to sample.
% nonempty_c_num is the number of non-empty clusters.
% n_c is of size [1, nonempty_c_num] and contains the number of patients in each cluster.
% t_SA is the coolness.

% Change Log:
%
%     1.1          30:aug:24    mx243      First version.
%     1.16         03:sep;24    rfs34      Parameter descriptions added and unnecessary min removed.
%     1.23         11:sep:24    rfs34      Fixed bug due to not accounting for reordering to put
%                                          non-empty clusters first always.
%     1.25         12:sep:24    mx243      Fixed bug due to forgetting that sum out of date.
%     1.29         16:sep:24    rfs34      Now assumes c annealed also.
%     1.30         16:sep:24    rfs34      Clips prior on C at max_C.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_largeC.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

lb_C_1 = max(floor(C / 2), nonempty_c_num); % The lb and ub of possible values of C' given C.
ub_C_1 = min(2 * C + 1, Priors.max_C);

log_Chi = NaN(1, 4 * C + 3);
log_norm_const_1 = -Inf; % Normalisation const for P' ( C' | C ).

for it_C = lb_C_1 : ub_C_1 % Compute log_Chi for P' ( C' | C ).
    sum = (it_C - 1) * log(Priors.kappa_C) - nonempty_c_num * gammaln(Priors.kappa_eta / it_C);
    for it_c = 1 : nonempty_c_num
        sum = sum + gammaln(Priors.kappa_eta / it_C + t_SA * n_c(it_c));
    end
    log_Chi(it_C) = sum + gammaln(it_C + 1) - gammaln(it_C - nonempty_c_num + 1);
    log_norm_const_1 = -Eadd(-log_norm_const_1, -log_Chi(it_C));
end

% Note sample_finite_discrete() only returns the index.
C_prime = sample_finite_discrete(log_Chi(lb_C_1 : ub_C_1)) + lb_C_1 - 1;

lb_C_2 = max(floor(C_prime / 2), nonempty_c_num); % The lb and ub of possible values of C given C'.
ub_C_2 = min(2 * C_prime + 1, Priors.max_C);

for it_C = lb_C_2 : ub_C_2 % Compute log_Chi for P' ( C | C' ).
    if it_C < lb_C_1 || it_C > ub_C_1 % Not computed before.
        sum = (it_C - 1) * log(Priors.kappa_C) - nonempty_c_num * gammaln(Priors.kappa_eta / it_C);
        for it_c = 1 : nonempty_c_num
            sum = sum + gammaln(Priors.kappa_eta / it_C + t_SA * n_c(it_c));
        end
        log_Chi(it_C) = sum + gammaln(it_C + 1) - gammaln(it_C - nonempty_c_num + 1);
    end
end

log_norm_const_2 = -Esum(-log_Chi(lb_C_2 : ub_C_2), 2); % Normalisation const for P' ( C | C' ).

% 2nd entry: log_Chi(C_prime) + log_Chi(C) - log_norm_const_2 - log_Chi(C)
% - log_C(C_prime) + log_norm_const_1 = log_norm_const_1 - log_norm_const_2.
logA = log_norm_const_1 - log_norm_const_2;
u = rand;
if log(u) < logA
    C = C_prime;
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
