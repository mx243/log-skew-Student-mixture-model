function [Samples, synth_seed] = get_prior_smpl(Priors, Nsmpl, synth_seed)
% Usage: [Samples, synth_seed] = get_prior_smpl(Priors, Nsmpl, synth_seed)
% This function returns Nsmpl naive prior samples from P(\theta).

% Change Log:
%
%     1.1          08:oct:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.get_prior_smpl.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 3
    synth_seed = [];
end

if ~isempty(synth_seed) % Save previous settings if given synth seed, then set rng using this seed.
    previousSettings = rng(synth_seed);
end

tmp_synth_seed = rng; % Record the settings used to synthesize this set of data.

D = Priors.D;
K = 0;

Samples = struct('x1', NaN(K, Nsmpl), 'alpha', NaN(K, Nsmpl), 'c', NaN(K, Nsmpl), 'nonempty_c_num', NaN(1, Nsmpl), ...
                     'C', NaN(1, Nsmpl), 'phi', struct('p', cell(1, Nsmpl), 'm', cell(1, Nsmpl), ...
                                                   'mu', cell(1, Nsmpl), 'S', cell(1, Nsmpl), 'nu', cell(1, Nsmpl), ...
                                                   'n_c', cell(1, Nsmpl), 'ks_at_c', cell(1, Nsmpl)), ...
                     'S_mu', NaN(D, D, Nsmpl), 'mu_mu', NaN(D, Nsmpl), 'm_S', NaN(1, Nsmpl), ...
                     'R_S', NaN(D, D, Nsmpl), 'R_S_mu', NaN(D, D, Nsmpl), ...
                     't_SA', NaN(1, Nsmpl), 'LN0', NaN, 'LN1', NaN, 'is_prior', NaN);

for it_N = 1 : Nsmpl
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

    % Record the sampled values.
    Samples.C(it_N) = C; 
    Samples.phi(it_N).p = p';
    Samples.phi(it_N).m = m';
    Samples.phi(it_N).mu = mu;
    Samples.phi(it_N).S = S;
    Samples.phi(it_N).nu = nu;
    Samples.S_mu(:, :, it_N) = S_mu;
    Samples.mu_mu(:, it_N) = mu_mu;
    Samples.m_S(it_N) = m_S;
    Samples.R_S(:, :, it_N) = R_S;
    Samples.R_S_mu(:, :, it_N) = R_S_mu;        
end

Samples.LN0 = 1;
Samples.LN1 = Nsmpl;
Samples.is_prior = 1; % Indicator of prior samples = 1.

if ~isempty(synth_seed) % Restore previous settings if rng was set using given synth seed.
     rng(previousSettings);
end

synth_seed = tmp_synth_seed; % Saved as the seed used in this run.

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
