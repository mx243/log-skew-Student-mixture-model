function smpl_file = run_get_prior_smpl(prior_file, Nsmpl, save_dir, save_id, ...
                                        mode, data_file, ctrl_file, synth_seed)
% Usage: smpl_file = run_get_prior_smpl(prior_file, Nsmpl, save_dir, save_id, ...
%                                       mode, data_file, ctrl_file, synth_seed)
% This function draws Nsample prior samples. When mode = 0, the samples are
% drawn from the naive prior P(\theta). When mode = 1, the samples are
% drawn from the prior P(\theta | all x_{ob \bar{1}}).

% Change Log:
%
%     1.1          08:oct:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.run_get_prior_smpl.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 5
    mode = 0;
end

if nargin < 8
    synth_seed = [];
end

if mode
    % Run sample_Stu_mix with x_ob in data_file and forcecensor = 1.
    smpl_file = sample_Stu_mix(data_file, prior_file, Nsmpl, save_dir, save_id, ctrl_file, [], synth_seed, 1);
    M = matfile(smpl_file, 'Writable', true);
    Samples = M.Samples;
    Samples.is_prior = 1; % Indicator of prior samples is set to 1.
    M.Samples = Samples;
else
    smpl_file = sprintf('%s/Stu_mix_smpls.%s.mat', save_dir, save_id);
    Priors = load(prior_file, 'kappa_C', 'max_C', 'kappa_eta', 'a_m', 'b_m', 'N_m', 'typ_m', 'kappa_nu', ...
                          'm_S_mu', 'm_R_S_mu', 'R_R_S_mu', 'mu_mu_mu', 'S_mu_mu', ...
                          'a_m_S', 'b_m_S', 'D', 'm_R_S', 'R_R_S');
    [Samples, synth_seed] = get_prior_smpl(Priors, Nsmpl, synth_seed);
    save(smpl_file, '-v6');
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
