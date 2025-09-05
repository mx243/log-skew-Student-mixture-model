function smpl_file = sample_Stu_mix(data_file, prior_file, smpl_num, save_dir, save_id, ...
                                    ctrl_file, Enables, analysis_seed, forcecensor, is_ASI)
% Usage: smpl_file = sample_Stu_mix(data_file, prior_file, smpl_num, save_dir, save_id, ...
%                                   ctrl_file, Enables, analysis_seed, forcecensor, is_ASI)
% This function draws smpl_num samples of \theta from its posterior in the
% log-skew-Student mixture model using MCMC and SA.
% The function also calculates the log of the evidence integral using 
% thermodynamic integration during the SA process, and return the value as
% 'I_thermo'.
% 'Samples' contains all the samples of \theta we drew in the entire
% MCMC plus annealing process. To extract the samples we want, access the
% last 'smpl_num' samples in 'Samples'.
% The input x_ob (D * K matrix) needs to be processed first: record the k's
%  st. x_ob(1, k) < 0 (censored patients), then set x_ob(1, k) = 
%  log(abs(x_ob(1, k))). 
% The *other* entries of x_ob are now expected to have had any transformations
%  such as logarithms or squares of logarithms done already.
% 
% data_file may be either a file containing x_ob (the observed data) and maybe Args_true, 
%  or an array which is itself the observed data.
%
% Below is the description for the Ctrl variables. This is a feature
% inherited from 'SA_MCMC_sample_Student.m' in Ex6 v1.13.
% Ctrl.t_SA_num: number of coolness values t_SA
% Ctrl.t_SA: a sequence of size t_SA_num containing coolness values t_SA
%            increasing from 0 to 1.
% Ctrl.nsmpl: a sequence of size t_SA_num indicating how many palindromic cycles
%             we do at each coolness.
% Ctrl.dis: a sequence of size t_SA_num indicating how many samples we discard at
%           each coolness when calculating the thermodynamic integral.
% Ctrl.dirsign: a sequence from {-1, 0, +1} indicating direction of movement of coolness.
% Ctrl.Nschedinsert: the coolness sequence index at which samples at t_SA = 1 can be added.
% The last (nsmpl(Nschedinsert) - dis(Nschedinsert)) samples drawn at t_SA = 1 are considered to 
% be examples from the target posterior distribution. 
% We will draw extra samples at t_SA = 1 if necessary, to make sure the
% last smpl_num samples are all drawn at t_SA = 1. These will be returned
% as Samples.
% N will denote the total number of palindromic cycles we did. All of these
% samples will be returned as Inter_smpls.
% If forcecensor is passed and is non-zero then all x_ob(1, :) will be forced to be censored at
%  -Inf, but only after S_0 and mu_0 have been determined. The purpose of this is to produce
%  a value to be subtracted from I_thermo to get P(x(1, :) | x(2 : end, :)).
% If is_ASI = 1, then we do not log x1. This is used for calculating
% ASI.

% Change Log:
%
%     1.1          22:aug:24    mx243      First version.
%     1.15         02:sep:24    mx243      log all entries of x and
%                                          preallocate space for the 
%                                          samples using arrays of cells.
%     1.16         03:sep:24    rfs34      No longer takes logs of all but the first entry of x_ob,
%                                          as some entries may need other transformations 
%                                          such as squaring after taking logarithms, so this 
%                                          needs to be done before calling this function;
%                                          vectorised initial calcs; now checks validity of Enables;
%                                          indexing of structures fixed; missing checks of Enables added;
%                                          c passed to resample_alpha.
%     1.17         03:sep:24    mx243      Added missing resampling steps for R_S_mu in MCMC.
%                                          Added description for the Ctrl
%                                          variables.
%     1.18         04:sep:24    mx243      Added nonempty_c_num to Samples, Inter_smpls.
%                                          Change undefined r to m on line 222.
%     1.19         06:sep:24    mx243      Changed the content of Samples,
%                                          incorporated the wrapper function 
%                                          run_Stu_mix, saves everything.
%     1.20         09:sep:24    rfs34      x1 added to Args_true; \r changed to \n at the end,
%                                          which matters on Linux.
%     1.21         10:sep:24    rfs34      Set the initial value of x1 to
%                                          be Args_true.x1 the uncensored true value.
%     1.22         10:sep:24    rfs34      Now allow censoring to stay on if x1 is disabled.
%     1.23         10:sep:24    mx243      Added return value smpl_file, added -v6 to the save command;
%                                          rfs34: also allowed censoring to remain if x1 disabled.
%     1.26         13:sep:24    rfs34      Now takes typ_m from Priors and uses it.
%     1.28         13:sep:24    rfs34      Now remembers to pass r to resample_c as well.
%     1.29         16:sep:24    rfs34      Now assumes c is being annealed also.
%     1.30         16:sep:24    rfs34      Clips prior on C at max_C.
%     1.32         16:sep:24    rfs34      forcecensor option added.
%     1.33         17:sep:24    rfs34      Changed method of calculating I_thermo when c annealed.
%     1.34         17:sep:24    rfs34      Added p to the parameter list of F_of_theta_and_x.
%     1.36         18:sep:24    rfs34      Fixed setting of mu_0 and S_0 when size(x,1) <= 2.
%     1.37         18:sep:24    mx243      Added is_ASI.
%     1.38         20:sep:24    rfs34      Trivial bug fixed.
%     1.41         21:sep:24    rfs34      Separates I_thermo in two components and increased sigma_0;
%                                          allows for annealing backwards.
%     1.43         24:sep:24    rfs34      Management of Ctrl updated; samples at t_SA=1 are now
%                                          in indices LN0 : LN1.
%     1.45         03:oct:24    rfs34      Merged with changes made by mx243; in particular 
%                               mx243      data_file can now be the array x_ob rather than a filename. 
%     1.50         08:oct:24    mx243      Added LN0, LN1, is_prior to Samples.
%     1.53         16:oct:24    rfs34      Now assumes is_synth=0 if not in the data file.
%     1.57         21:oct:24    mx243      Added a resampling step which
%                                          moves mu1 and x1 together in order to increase mobility.
%     1.59         25:oct:24    mx243      Added coefficient for the correct evidence integral.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_Stu_mix.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

Priors = load(prior_file, 'kappa_C', 'max_C', 'kappa_eta', 'a_m', 'b_m', 'N_m', 'typ_m', 'kappa_nu', ...
                          'm_S_mu', 'm_R_S_mu', 'R_R_S_mu', 'mu_mu_mu', 'S_mu_mu', ...
                          'a_m_S', 'b_m_S', 'D', 'm_R_S', 'R_R_S');
D = Priors.D;

if nargin < 6 || isempty(ctrl_file)
    Ctrl = struct('t_SA_num', 1, 't_SA', 1, 'nsmpl', smpl_num, 'dis', 0, ...
                  'dirsign', 0, 'Nschedinsert', 1); % No annealing.
else
    Ctrl = load(ctrl_file, 't_SA_num', 't_SA', 'nsmpl', 'dis', 'dirsign', 'Nschedinsert');
    if(isempty(Ctrl.nsmpl)) % Default annealing: draw one sample at each t_SA, discard none.
        Ctrl.nsmpl = ones(1, Ctrl.t_SA_num);
        Ctrl.dis = zeros(1, Ctrl.t_SA_num);
    end
end

if nargin < 7 || isempty(Enables),
    Enables = struct('x1', 2, 'alpha', 2, 'c', 2, 'phi', struct('C', 2, 'p', 2, 'm', 2, 'mu', 2, 'S', 2, 'nu', 2), ...
                     'S_mu', 2, 'mu_mu', 2, 'm_S', 2, 'R_S', 2, 'R_S_mu', 2);
end

if nargin < 8,
    analysis_seed = [];
end

if nargin < 9 || isempty(forcecensor),
    forcecensor = 0;
end

if nargin < 10 || isempty(is_ASI),
    is_ASI = 0;
end

if ~isempty(analysis_seed) % Save previous settings if given analysis seed, then set rng using this seed.
    previousSettings = rng(analysis_seed);
end

tmp_analysis_seed = rng; % Record the settings of rng before the run starts.

if ischar(data_file)
    load(data_file, 'x_ob', 'is_synth');
    if ~exist('is_synth', 'var'),
       is_synth = 0;
    end
else
    x_ob = data_file; % If data_file is not a char vector, the it's the data itself.
    is_synth = 0; % Does not take in Args_true in this case even if the data is synthetic.
end

K = size(x_ob, 2);
if is_synth
    load(data_file, 'Args_true');
else
    % Value given to phi here doesn't matter as long as it can pick up errors.
    Args_true = struct('x1', NaN(1, K), 'alpha', NaN(1, K), 'c', NaN(1, K), ...
                       'phi', struct('C', 0, 'p', zeros(1, 0), 'm', zeros(1, 0), ...
                                     'mu', zeros(D, 0), 'S', zeros(D, D, 0), 'nu', zeros(D, 0)), ...
                       'S_mu', NaN(D, D), 'mu_mu', NaN(D, 1), 'm_S', NaN, 'R_S', NaN(D, D), 'R_S_mu', NaN(D, D));
end

if size(x_ob, 1) ~= D,
    error('sample_Stu_mix called with Prior with wrong dimension.')
end

% Get parameters for the proposed distribution P_0 of log(x1).
% They are now based on both deaths and censoring times,
% where in previous versions they were based only on censoring times.
x1 = x_ob(1, :);

if isempty(x1),
   mu_0 = 0;
   sigma_0 = 1;
elseif isscalar(x1),
   if is_ASI,
      mu_0 = x1;
   else
      mu_0 = log(abs(x1));
   end
   sigma_0 = 1;
else
   if is_ASI,
      mu_0 = mean(x1);
      sigma_0 = std(x1);
   else
      mu_0 = mean(log(abs(x1)));
      sigma_0 = std(log(abs(x1)));
   end
end
if sigma_0 == 0,
    warning('All death and censoring times are the same; sigma_0 is 0 and I_thermo will be wrong.');
end

sigma_0 = 5 * sigma_0; 

S_0 = 1 / (sigma_0 ^ 2);

% Process x1.
if is_ASI
    is_censored = false(1, K); % Used to calculate ASI, no censoring.
else
    is_censored = x_ob(1, :) < 0;
    if forcecensor,
        is_censored(:) = true;
    end
    x_ob(1, :) = log(abs(x_ob(1, :)));
    if forcecensor,
        x_ob(1, :) = -Inf;
    end
end
logd = x_ob(1, :);

% Sum of logd of uncensored patients.
evidence_coeff = sum(logd(~is_censored), 2);

% Total number of samples we need to draw.
Nextraatone = max(0, smpl_num - (Ctrl.nsmpl(Ctrl.Nschedinsert) - Ctrl.dis(Ctrl.Nschedinsert)));
Ctrl.nsmpl(Ctrl.Nschedinsert) = Ctrl.nsmpl(Ctrl.Nschedinsert) + Nextraatone;
N = sum(Ctrl.nsmpl, 2); % Total number of samples to be recorded.

% Moved C out from phi for easier access and added supplimentary n_c and
% ks_at_c (not actually part of phi, just for convenience).
% Also since we're now saving smpl_num too, simply extract the last
% smpl_num samples in Samples to get the samples we want (at t_SA = 1).
Samples = struct('x1', NaN(K, N), 'alpha', NaN(K, N), 'c', NaN(K, N), 'nonempty_c_num', NaN(1, N), ...
                     'C', NaN(1, N), 'phi', struct('p', cell(1, N), 'm', cell(1, N), ...
                                                   'mu', cell(1, N), 'S', cell(1, N), 'nu', cell(1, N), ...
                                                   'n_c', cell(1, N), 'ks_at_c', cell(1, N)), ...
                     'S_mu', NaN(D, D, N), 'mu_mu', NaN(D, N), 'm_S', NaN(1, N), ...
                     'R_S', NaN(D, D, N), 'R_S_mu', NaN(D, D, N), ...
                     't_SA', NaN(1, N), 'LN0', NaN, 'LN1', NaN, 'is_prior', NaN);

% Check validity of Enables.
if Enables.phi.C == 2,
   if ~all([Enables.phi.p; Enables.phi.m; Enables.phi.mu; Enables.phi.S; Enables.phi.nu] == 2),
      error('Must not enable phi.C to 2 without enabling the rest of phi to 2 as well');
   end
end
if Enables.phi.C == 1,
   if ~all([Enables.phi.p; Enables.phi.m; Enables.phi.mu; Enables.phi.S; Enables.phi.nu] >= 1),
      error('Must not enable phi.C without enabling the rest of phi as well');
   end
end

% Assign initial values.
if Enables.R_S_mu == 2,
    R_S_mu = sample_Wis(Priors.m_R_S_mu, Priors.R_R_S_mu, D);
else
    R_S_mu = Args_true.R_S_mu;
end

if Enables.S_mu == 2,
    S_mu = sample_Wis(Priors.m_S_mu, R_S_mu, D);
else
    S_mu = Args_true.S_mu;
end

if Enables.mu_mu == 2,
    mu_mu = sample_Normal(Priors.mu_mu_mu, Priors.S_mu_mu);
else
    mu_mu = Args_true.mu_mu;
end

if Enables.m_S == 2,
    m_S = sample_proGamma(Priors.a_m_S, Priors.b_m_S, D, 1);
else
    m_S = Args_true.m_S;
end

if Enables.R_S == 2,
    R_S = sample_Wis(Priors.m_R_S, Priors.R_R_S, D);
else
    R_S = Args_true.R_S;
end

if Enables.phi.C == 2,
    if Enables.c == 2
        C = sample_discrete_Exp(Priors.kappa_C, 1, Priors.max_C);
    else % If c starts from a fixed value then must have C >= max_k{c_k}.
        C = sample_discrete_Exp(Priors.kappa_C, 1, Priors.max_C - Args_true.nonempty_c_num) ...
            + Args_true.nonempty_c_num - 1;
    end
else
    C = Args_true.phi.C;
end
eta = (Priors.kappa_eta / C) * ones(1, C);

if Enables.phi.p == 2
    p = sample_Dir(eta, C);
else
    p = Args_true.phi.p;
end

if Enables.phi.m == 2
    m = NaN(1, C);
    for c = 1 : C % Have to do this because ars.m can't draw more than one samples at once.
        m(c) = sample_proGamma(Priors.a_m, Priors.b_m, Priors.N_m, Priors.typ_m); % In fact N_m = 1.
    end
else
    m = Args_true.phi.m;
end
r = r_from_m(m, Priors.typ_m);

if Enables.phi.mu == 2
    mu = sample_Normal(mu_mu, S_mu, C);
else
    mu = Args_true.phi.mu;
end

if Enables.phi.S == 2
    S = NaN(D, D, C);
    for c = 1 : C 
        S(:, :, c) = sample_Wis(m_S, (m_S - 1) * R_S, D);
    end
else
    S = Args_true.phi.S;
end

if Enables.phi.nu == 2
    nu = NaN(D, C);
    for c = 1 : C 
        nu(:, c) = sample_Normal(zeros(D, 1), S(:, :, c) / Priors.kappa_nu);
    end
else
    nu = Args_true.phi.nu;
end

if Enables.c == 2 
    c = sample_finite_discrete(log(p), K);
else
    c = Args_true.c;
end
[c, n_c, ks_at_c, nonempty_c_num, p, m, mu, S, nu] = reorder_col(C, K, c, p, m, mu, S, nu);

if Enables.alpha == 2
    alpha = NaN(1, K);
    for k = 1 : K
        alpha(k) = sample_Gamma(m(c(k)), r_from_m(m(c(k)), Priors.typ_m));
    end
else
    alpha = Args_true.alpha;
end

x = x_ob;
if Enables.x1 == 2,
    for k = 1 : K
        if is_censored(k)
            x(1, k) = sample_truncated_Normal(mu_0, S_0, logd(k)); % Distr of x1 when t_SA = 0.
        end
    end
else
   x(1, :) = Args_true.x1;
   % We allow censoring to remain in place here, as x1 is (perhaps partially) disabled,
   %  but x(1, :) has its correct true value (perhaps only to start with) anyhow.
end

% End of initialisation.

smpl_cnt = 0;
E_F_of_theta_and_x = NaN(2, Ctrl.t_SA_num); % The value of the integrand at each t_SA.

xmuaccepts = 0;
xmutries = 0;

% The MCMC process.
for j = 1 : Ctrl.t_SA_num
    t_SA = Ctrl.t_SA(j);
    tmp_E_F_of_theta_and_x = [0; 0];

    for it = 1 : Ctrl.nsmpl(j)

        % Resample x1.
        if(Enables.x1)
            x = Stu_mix_resample_x1(x, alpha, D, K, mu, S, nu, c, logd, is_censored, mu_0, S_0, t_SA);
        end

        % Resample alpha.
        if(Enables.alpha)
            alpha = Stu_mix_resample_alpha(alpha, x, K, mu, S, nu, m, r, c, t_SA);
        end

        % Resample c.
        if(Enables.c) 
            c = Stu_mix_resample_c(c, x, alpha, K, p, mu, S, nu, m, r, C, t_SA);
            [c, n_c, ks_at_c, nonempty_c_num, p, m, mu, S, nu] = reorder_col(C, K, c, p, m, mu, S, nu);
        end

        % Resample phi.
        [C, p, m, mu, S, nu] = Stu_mix_resample_phi(C, p, m, mu, S, nu, ...
                                                    x, alpha, K, c, n_c, ks_at_c, nonempty_c_num, ...
                                                    Priors, S_mu, mu_mu, m_S, R_S, D, t_SA, Enables);
        r = r_from_m(m, Priors.typ_m);

        % Extra resampling step that moves x1 and mu1 together.
        if(Enables.x1 && Enables.phi.mu)
            [x, mu, xmuaccepts, xmutries] ...
               = Stu_mix_resample_x1_with_mu1(x, mu, D, nonempty_c_num, n_c, ks_at_c, ...
                                              is_censored, logd,  mu_0, S_0, mu_mu, S_mu, t_SA, xmuaccepts, xmutries);
        end
        
        % Resample mu_mu.
        if (Enables.mu_mu),
           mu_mu = Stu_mix_resample_mu_mu(Priors, mu, C, S_mu);
        end

        % Resample S_mu.
        if (Enables.S_mu),
           S_mu = Stu_mix_resample_S_mu(Priors, mu, D, C, mu_mu, R_S_mu);
        end

        % Resample R_S_mu.
        if (Enables.R_S_mu),
           R_S_mu = Stu_mix_resample_R_S_mu(Priors, D, S_mu);
        end

        % Resample m_S.
        if (Enables.m_S),
           m_S = Stu_mix_resample_m_S(Priors, D, C, S, R_S);
        end

        % Resample R_S. This is the middle of the palindromic cycle.
        if Enables.R_S,
           R_S = Stu_mix_resample_R_S(Priors, D, C, S, m_S);
        end

        % Resample m_S.
        if Enables.m_S,
           m_S = Stu_mix_resample_m_S(Priors, D, C, S, R_S);
        end

        % Resample R_S_mu.
        if (Enables.R_S_mu),
           R_S_mu = Stu_mix_resample_R_S_mu(Priors, D, S_mu);
        end

        % Resample S_mu.
        if Enables.S_mu,
           S_mu = Stu_mix_resample_S_mu(Priors, mu, D, C, mu_mu, R_S_mu);
        end

        % Resample mu_mu.
        if Enables.mu_mu,
           mu_mu = Stu_mix_resample_mu_mu(Priors, mu, C, S_mu);
        end

        % Extra resampling step that moves x1 and mu1 together.
        if(Enables.x1 && Enables.phi.mu)
            [x, mu, xmuaccepts, xmutries] ...
               = Stu_mix_resample_x1_with_mu1(x, mu, D, nonempty_c_num, n_c, ks_at_c, ...
                                              is_censored, logd,  mu_0, S_0, mu_mu, S_mu, t_SA, xmuaccepts, xmutries);
        end

        if any(x(1, :) < logd),
            error('Constraint violated.');
        end

        % Resample phi.
        [C, p, m, mu, S, nu] = Stu_mix_resample_phi(C, p, m, mu, S, nu, ...
                                                    x, alpha, K, c, n_c, ks_at_c, nonempty_c_num, ...
                                                    Priors, S_mu, mu_mu, m_S, R_S, D, t_SA, Enables);
        r = r_from_m(m, Priors.typ_m);

        % Resample c.
        if(Enables.c) 
            c = Stu_mix_resample_c(c, x, alpha, K, p, mu, S, nu, m, r, C, t_SA);
            [c, n_c, ks_at_c, nonempty_c_num, p, m, mu, S, nu] = reorder_col(C, K, c, p, m, mu, S, nu);
        end

        % Resample alpha.
        if(Enables.alpha)
            alpha = Stu_mix_resample_alpha(alpha, x, K, mu, S, nu, m, r, c, t_SA);
        end

        % Resample x1.
        if(Enables.x1)
            x = Stu_mix_resample_x1(x, alpha, D, K, mu, S, nu, c, logd, is_censored, mu_0, S_0, t_SA);
        end

        % End of the palindromic cycle.

        % Samples drawn when it <= Ctrl.dis(j) are ignored in the
        % calculation of E(F(theta, x)).
        if(it > Ctrl.dis(j))
            tmp_E_F_of_theta_and_x = tmp_E_F_of_theta_and_x + ...
                                     F_of_theta_and_x(x, alpha, c, p, nonempty_c_num, ...
                                                      mu, S, nu, is_censored, logd, mu_0, S_0, D, K);
        end

        smpl_cnt = smpl_cnt + 1;

        % Record the resampled values.
        Samples.x1(:, smpl_cnt) = x(1, :)';
        Samples.alpha(:, smpl_cnt) = alpha';
        Samples.c(:, smpl_cnt) = c';
        Samples.nonempty_c_num(smpl_cnt) = nonempty_c_num;
        Samples.C(smpl_cnt) = C;

        % Because of how Inter_smpls.phi has been defined, it is an array of structures
        % not a structure of arrays.
        Samples.phi(smpl_cnt).p = p';
        Samples.phi(smpl_cnt).m = m';
        Samples.phi(smpl_cnt).mu = mu;
        Samples.phi(smpl_cnt).S = S;
        Samples.phi(smpl_cnt).nu = nu;
        Samples.phi(smpl_cnt).n_c = n_c;
        Samples.phi(smpl_cnt).ks_at_c = ks_at_c;

        Samples.S_mu(:, :, smpl_cnt) = S_mu;
        Samples.mu_mu(:, smpl_cnt) = mu_mu;
        Samples.m_S(smpl_cnt) = m_S;
        Samples.R_S(:, :, smpl_cnt) = R_S;
        Samples.R_S_mu(:, :, smpl_cnt) = R_S_mu;
        Samples.t_SA(smpl_cnt) = t_SA;

        fprintf('\rMCMC: ...done %d of %d, x1 with mu1 acceptrate = %g (%d of %d) ...', ...
                smpl_cnt, N, (xmuaccepts + 1) / (xmutries + 2), xmuaccepts, xmutries);
    end
    % Value of the integrand at t_SA.
    E_F_of_theta_and_x(:, j) = tmp_E_F_of_theta_and_x / (Ctrl.nsmpl(j) - Ctrl.dis(j));
end
fprintf('\n');

% Numerical integration.
if(Ctrl.t_SA_num > 1)
    I_thermo = trapz(Ctrl.t_SA, Ctrl.dirsign([1; 1], :) .* E_F_of_theta_and_x, 2);
    I_thermo(1, 1) = I_thermo(1, 1) - evidence_coeff;
else
    I_thermo = [0; 0];
    warning('I_thermo is calculated without annealing.');
end

LN0 = sum(Ctrl.nsmpl(1 : Ctrl.Nschedinsert - 1), 2) + Ctrl.dis(Ctrl.Nschedinsert) + 1;
LN1 = LN0 + smpl_num - 1;

Samples.LN0 = LN0;
Samples.LN1 = LN1;
Samples.is_prior = 0; % Not prior samples.

% For easier access. Access LN0 to LN1 in Samples (or increase LN0 to discard potentially non-converged samples).

if ~isempty(analysis_seed) % Restore previous settings if rng was set using given analysis seed.
     rng(previousSettings);
end

analysis_seed = tmp_analysis_seed; % Saved as the seed used in this run.

if ~is_ASI
    smpl_file = sprintf('%s/Stu_mix_smpls.%s.mat', save_dir, save_id);
else
    smpl_file = sprintf('%s/ASI_smpls.%s.mat', save_dir, save_id);
end
save(smpl_file, '-v7.3'); 

fprintf('Remember that I_thermo is log(P(x)) not log(P(x(1, :) | x(2 : end, :))),\n');
fprintf('  but the latter is what you want when comparing with tow24''s model;\n');
fprintf('Repeat run with all x(1, :) censored at -Inf (using forcecensor = 1) to get log(P(x(2 : end, :))).\n');

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
