function run_Stu_mix(data_file, prior_file, save_dir, save_id, ctrl_file, Enables, analysis_seed)
% Usage: run_Stu_mix(data_file, prior_file, save_dir, save_id, ctrl_file, Enables, analysis_seed)
% This function reads synthetic data from data_file, calls 'sample_Stu_mix'
% with prior_file and ctrl_file to get [Samples, Inter_smpls, I_thermo] and
% save these results at save_dir with save_id.
% It can reproduce a set of results when given analysis_seed.

% Change Log:
%
%     1.1          04:sep:24    mx243      First version.
%     1.20         09:sep:24    rfs34      Warning added that this function is now redundant.
%     1.30         16:sep:24    rfs34      Clips prior on C at max_C.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.run_Stu_mix.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

warning('This function is now redundant - use sample_Stu_mix.m instead');

if nargin < 6 || isempty(Enables),
    Enables = struct('x1', 2, 'alpha', 2, 'c', 2, 'phi', struct('C', 2, 'p', 2, 'm', 2, 'mu', 2, 'S', 2, 'nu', 2), ...
                     'S_mu', 2, 'mu_mu', 2, 'm_S', 2, 'R_S', 2, 'R_S_mu', 2);
end

if nargin < 7,
    analysis_seed = [];
end

if ~isempty(analysis_seed) % Save previous settings if given analysis seed, then set rng using this seed.
    previousSettings = rng(analysis_seed);
end

tmp_analysis_seed = rng; % Record the settings of rng before the run starts.

load(data_file, 'x_ob', 'Args_true'); 
Priors = load(prior_file, 'kappa_C', 'max_C', 'kappa_eta', 'a_m', 'b_m', 'N_m', 'typ_m', 'kappa_nu', ...
                          'm_S_mu', 'm_R_S_mu', 'R_R_S_mu', 'mu_mu_mu', 'S_mu_mu', ...
                          'a_m_S', 'b_m_S', 'D', 'm_R_S', 'R_R_S');
Ctrl = load(ctrl_file, 't_SA_num', 't_SA', 'nsmpl', 'dis', 'dirsign', 'Nschedinsert');
D = Priors.D;

% The aim of this function is to test convergence, so set smpl_num = 1.
[Samples, Inter_smpls, I_thermo] = sample_Stu_mix(x_ob, Priors, 1, Ctrl, Args_true, Enables);

if ~isempty(analysis_seed) % Restore previous settings if rng was set using given analysis seed.
     rng(previousSettings);
end

analysis_seed = tmp_analysis_seed; % Saved as the seed used in this run.

save(sprintf('%s/Stu_mix_smpls.%s.mat', save_dir, save_id)); 

return;

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
