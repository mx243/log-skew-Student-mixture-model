function prior_file = write_Stu_mix_priors(save_dir, save_id, mode, x)
% Usage: prior_file = write_Stu_mix_priors(save_dir, save_id, mode, x)
% This function creates a .mat file specifying a set of prior arguments for
% the Student mixture model.
% mode = -1: prior for modelling demeaned, descaled log(cytokine) values in the Whitworth dataset.
% mode = 0 : prior for modelling distribution of samples j of ASI.
% mode = 1 : prior for primary data analysis, with parameters
%            for top level of hierarchy taken from data x.
% mode = 2 : prior for analysis of Shariat prostate cancer dataset.

% Change Log:
%
%     1.1          04:sep:24    mx243      First version.
%     1.19         08:sep:24    mx243      Now can switch between prior for
%                                          ASI and for data analysis.
%     1.20         09:sep:24    rfs34      Comments changed to be clearer.
%     1.22         10:sep:24    mx243      Added return value prior_file, added -v6 to the save command.
%     1.26         13:sep:24    rfs34      Added typ_m to the priors.
%     1.30         16:sep:24    rfs34      Added max_C to the priors.
%     1.37         18:sep:24    mx243      D = size(x, 1) on line 56.
%     1.53         15:oct:24    mx243      Process x before using it to get the top level parameters.
%     1.54         17:oct:24    mx243      New set of parameters for mode 1, adapted to old dataset.
%     1.70         17:jan:25    rfs34      Added prior as for old Shariat dataset.
%     1.76         21:mar:25    mx243      Added mode for modelling cytokine values in the Whitworth dataset.
%     1.77         22:mar:25    rfs34      Adjusted cytokine priors to allow some narrower peaks.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.write_Stu_mix_priors.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

% Both cases are subject to change.
if nargin < 4
    if mode > 0
        error('write_Stu_mix_priors called with mode 1 or 2 and no data.');
    end
end

if mode == -1
    D = 1;
    N_m = 1;
    kappa_nu = 1;
    kappa_eta = 10;
    kappa_C = 0.85;
    max_C = 5;
    mu_mu_mu = 0;
    S_mu_mu = 1;
    m_S_mu = 1.1 / 2;
    m_R_S_mu = 2;
    R_R_S_mu = 2.8 / 1.5;
    a_m = 1;
    b_m = 3;
    typ_m = 1;
    m_R_S = 2 * 4;
    R_R_S = 200 / 16;
    a_m_S = 1;
    b_m_S = 2;
end

if mode == 0
    D = 1;
    N_m = 1;
    kappa_nu = 1;
    kappa_eta = 10;
    kappa_C = 0.9;
    max_C = 20;
    mu_mu_mu = 0;
    S_mu_mu = 1;
    m_S_mu = 1.1;
    m_R_S_mu = 2;
    R_R_S_mu = 2.8;
    a_m = 1;
    b_m = 3;
    typ_m = 1;
    m_R_S = 2;
    R_R_S = 200;
    a_m_S = 1;
    b_m_S = 2;

elseif mode == 1,

    % Process input x by applying log(abs()) to the first row.
    x(1, :) = log(abs(x(1, :)));

    D = size(x, 1);
    N_m = 1;
    kappa_nu = 0.256;
    kappa_eta = 10;
    kappa_C = 0.7;
    max_C = 7;
    mu_mu_mu = mean(x, 2);

    T = cov(x', 1) ^ (-1);

    S_mu_mu = 0.32 * T;
    m_S_mu = 2.5;
    m_R_S_mu = 2;
    R_R_S_mu = 12 * T;
    a_m = 1;
    b_m = 2;
    typ_m = 2;
    m_R_S = 71;
    R_R_S = 800 * T; 
    a_m_S = 1;
    b_m_S = 0.01;

elseif mode == 2,

    % Process input x by applying log(abs()) to the first row.
    x(1, :) = log(abs(x(1, :)));

    D = size(x, 1);
    N_m = 1;
    kappa_nu = 0.1;
    kappa_eta = 10;
    kappa_C = 0.09;
    max_C = 7;
    mu_mu_mu = mean(x, 2);

    T = cov(x', 1) ^ (-1);

    S_mu_mu = 6.25 * T;
    m_S_mu = 2.5;
    m_R_S_mu = 2;
    R_R_S_mu = 24 * T;
    a_m = 1;
    b_m = 20;
    typ_m = 2;
    m_R_S = 133;
    R_R_S = 140 * T; 
    a_m_S = 1;
    b_m_S = 0.01;

end

prior_file = sprintf('%s/priors_Stu_mix.%s.mat', save_dir, save_id);
save(sprintf('%s/priors_Stu_mix.%s.mat', save_dir, save_id), '-v6');

return;

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
