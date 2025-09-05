function [S, nu] = Stu_mix_resample_S_and_nu_empty(Priors, S, nu, m_S, R_S, nonempty_c_num, C, D, Enables)
% Usage: [S, nu] = Stu_mix_resample_S_and_nu_empty(Priors, S, nu, m_S, R_S, nonempty_c_num, C, D, Enables)
% This function resamples S, nu together at c's that are empty.

% Change Log:
%
%     1.1          05:sep:24    mx243      First version.
%     1.19         06:sep:24    mx243      Comment only changed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.Stu_mix_resample_S_and_nu_empty.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

% First S.
if Enables.phi.S
    for it_c = nonempty_c_num + 1 : C
        S(:, :, it_c) = sample_Wis(m_S, (m_S - 1) * R_S, D);
    end
end

% Then nu | S.
if Enables.phi.nu % In fact can't have C axis on and nu off. It is assumed Enables satisfies this requirement.
    for it_c = nonempty_c_num + 1 : C
        nu(:, it_c) = sample_Normal(zeros(D, 1), S(:, :, it_c) / Priors.kappa_nu); 
    end
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
