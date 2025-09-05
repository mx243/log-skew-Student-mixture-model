function m = Stu_mix_resample_m(Priors, alpha, n_c, ks_at_c, nonempty_c_num, C)
% Usage: m = Stu_mix_resample_m(Priors, alpha, n_c, ks_at_c, nonempty_c_num, C)
% This function resamples phi.mu in the Stu mixture model.
% nonempty_c_num is the number of non-empty clusters.
% n_c is of size [1, nonempty_c_num] and contains the number of patients in each cluster.
% ks_at_c has size [K, nonempty_c_num], where K is the number of patients;
%  it contains the list of patients in each cluster.

% Change Log:
%
%     1.1          31:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Comments only changed.
%     1.18         06:sep:24    mx243      Added correct final arguments to call to sample_proGamma.m .

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_m.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

m = NaN(1, C);
typ = Priors.typ_m;

for it_c = 1 : C
    if it_c <= nonempty_c_num

        tmp_sum = sum(alpha(ks_at_c(1 : n_c(it_c), it_c)'), 2) - ...
                  sum(log(alpha(ks_at_c(1 : n_c(it_c), it_c)')), 2) - n_c(it_c);
        m(it_c) = sample_proGamma(Priors.a_m + tmp_sum, Priors.b_m + n_c(it_c), 1, typ);
    else
        m(it_c) = sample_proGamma(Priors.a_m, Priors.b_m, 1, typ);
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
