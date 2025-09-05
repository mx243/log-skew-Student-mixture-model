function p = Stu_mix_resample_p(Priors, C, n_c, nonempty_c_num, t_SA)
% Usage: p = Stu_mix_resample_p(Priors, C, n_c, nonempty_c_num, t_SA)
% This function resamples p with given C in the Stu mixture model.
% nonempty_c_num is the number of non-empty clusters.
% n_c is of size [1, nonempty_c_num] and contains the number of patients in each cluster.

% Change Log:
%
%     1.1          30:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Comments only changed.
%     1.29         16:sep:24    rfs34      Now assumes c is being annealed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_p.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

lambda = (Priors.kappa_eta / C) * ones(1, C);
lambda(1 : nonempty_c_num) = lambda(1 : nonempty_c_num) + t_SA * n_c;
p = sample_Dir(lambda, C);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
