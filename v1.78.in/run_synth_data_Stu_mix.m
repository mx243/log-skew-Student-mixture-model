function data_file = run_synth_data_Stu_mix(prior_file, save_dir, save_id, K, mu_1, S_1, synth_seed)
% Usage: data_file = run_synth_data_Stu_mix(prior_file, save_dir, save_id, K, mu_1, S_1, synth_seed)
% This function calls synth_data_Stu_mix with prior parameters found at
% 'prior_file'. It saves the synthesized data as
% 'save_dir\data_Stu_mix_synth.save_id.mat'.
% If synth_seed is given, we use it to repeat the synthesis of a set of
% data.
% mu_1 and S_1 are the parameters for the Gaussian censoring distribution.

% Change Log:
%
%     1.1          04:sep:24    mx243      First version.
%     1.19         06:sep:24    mx243      Added is_synth. This indicator
%                                          will also be added to real data.
%     1.20         09:sep:24    rfs34      Comments only changed.
%     1.22         10:sep:24    mx243      Added return value data_file, added -v6 to the save command.
%     1.26         13:sep:24    rfs34      Now takes account of typ_m.
%     1.30         16:sep:24    rfs34      Clips prior on C at max_C.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.run_synth_data_Stu_mix.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 7
    synth_seed = [];
end

if ~isempty(synth_seed) % Save previous settings if given synth seed, then set rng using this seed.
    previousSettings = rng(synth_seed);
end

tmp_synth_seed = rng; % Record the settings used to synthesize this set of data.

Priors = load(prior_file, 'kappa_C', 'max_C', 'kappa_eta', 'a_m', 'b_m', 'N_m', 'typ_m', 'kappa_nu', ...
                          'm_S_mu', 'm_R_S_mu', 'R_R_S_mu', 'mu_mu_mu', 'S_mu_mu', ...
                          'a_m_S', 'b_m_S', 'D', 'm_R_S', 'R_R_S');

[x_ob, Args_true] = synth_data_Stu_mix(Priors, K, mu_1, S_1);
is_synth = 1;

if ~isempty(synth_seed) % Restore previous settings if rng was set using given synth seed.
     rng(previousSettings);
end

synth_seed = tmp_synth_seed; % Saved as the seed used in this run.

data_file = sprintf('%s/data_Stu_mix_synth.%s.mat', save_dir, save_id);
save(sprintf('%s/data_Stu_mix_synth.%s.mat', save_dir, save_id), '-v6');

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
