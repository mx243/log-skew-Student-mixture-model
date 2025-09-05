function m_S = Stu_mix_resample_m_S(Priors, D, C, S, R_S)
% Usage: m_S = Stu_mix_resample_m_S(Priors, D, C, S, R_S)
% This function resamples m_S in the Stu mixture model.

% Change Log:
%
%     1.1          1:sept:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_m_S.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

tmp_sum = trace(R_S * sum(S, 3)) + C * Edet(R_S) - D * C; % Edet(A) = -log(det(A)).
for it_c = 1 : C
    tmp_sum = tmp_sum + Edet(S(:, :, it_c));
end
m_S = sample_proGamma(Priors.a_m_S + tmp_sum, Priors.b_m_S + C, D, 1);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
