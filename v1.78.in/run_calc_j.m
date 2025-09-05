function ASI_data_file = run_calc_j(data_file, primary_prior_file, primary_smpl_num, primary_ctrl_file, ...
                                    save_dir, primary_id, ASI_id, propunseen, Nsplit, analysis_seed, Nskip)
% Usage: ASI_data_file = run_calc_j(data_file, primary_prior_file, primary_smpl_num, primary_ctrl_file, ...
%                                   save_dir, primary_id, ASI_id, propunseen, Nsplit, analysis_seed, Nskip)
% This function separates the data x into training and unseen data Nsplit 
% times, uses the tranining data to train a model and uses the model to get
% a vector of j on the unseen data. Then all the values of j are
% combined into one vector and saved in ASI_data_file.
% propunseen is the fraction of the patients that should be deemed unseen at each of Nsplit splits.

% Change Log:
%
%     1.1          29:sep:24    mx243      First version.
%     1.45         03:oct:24    rfs34      Layout only changed.
%     1.46         03:oct:24    rfs34      Made propunseen and Nsplit parameters, and changed 
%                                          filenames to be more useful.
%     1.50         07:oct:24    mx243      Trivial error fixed; each set of unseen data saved separately.
%     1.62         28:oct:24    mx243      Added Nskip option.
%     1.64         30:oct:24    rfs34      Added ability to resume from tempsave file.
%     1.65         31:oct:24    rfs34      Changed use of analysis_seed to cover splitting also.
%     1.73         26:jan:25    rfs34      Now reuses existing modelling runs and unseen sets
%                                          if available.
%     1.74         26:jan:25    rfs34      Minor bugs in above fixed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.run_calc_j.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 8 || isempty(propunseen),
   propunseen = exp(-1);
end
if nargin < 9 || isempty(Nsplit),
   Nsplit = 2;
end
if nargin < 10 || isempty(analysis_seed),
    analysis_seed = [];
end
if nargin < 11 || isempty(Nskip),
    Nskip = 1;
end

if ~isempty(analysis_seed) % Save previous settings if given analysis seed, then set rng using this seed.
    previousSettings = rng(analysis_seed);
end

tmp_analysis_seed = rng; % Record the settings of rng before the run starts.

load(data_file, 'x_ob');
K = size(x_ob, 2);
D = size(x_ob, 1);

N_unseen = ceil(propunseen * K);
j_vec = NaN(1, Nsplit * N_unseen);

% Infer the mean posterior lambda from the set of *all* death and censoring times,
% in order to make sure we are using the same value as others using different splits.

% Gamma prior on lambda, now set to match tow24's values.
m2 = 0.1;
r2 = 40; % days (so prior mean is 400 days)

% Posterior of lambda ~ Gamma(m2 + n_uncensored, r2 + sum(|x_ob(1, :)|)).
% We extract its mean which is then used for P2 in calc_j.m .
n_uncensored = sum(x_ob(1, :) > 0);
lambda = (m2 + n_uncensored) / (r2 + sum(abs(x_ob(1, :))));

it_N = 0;

if exist(sprintf('%s/tempsave.%s.mat', save_dir, primary_id), 'file'),
   load(sprintf('%s/tempsave.%s.mat', save_dir, primary_id));
   rng(tempsaveseed);
   clear tempsaveseed
end

% Get all values of log(j);
while it_N < Nsplit,
    it_N = it_N + 1;
 
    is_unseen = false(1, K);
    ind = randperm(K);
    is_unseen(ind(1 : N_unseen)) = true;

    training_data = x_ob(:, ~is_unseen); % Not passed to file.
    save_id = [primary_id, sprintf('.%d', it_N)];

    if (exist(sprintf('%s/Stu_mix_smpls.%s.mat', save_dir, save_id), 'file') ...
        && exist(sprintf('%s/x_ob_unseen.%s.mat', save_dir, save_id), 'file')),

        warning(sprintf('Reusing old files for it_N = %d', it_N));
        trained_smpl_file = sprintf('%s/Stu_mix_smpls.%s.mat', save_dir, save_id);
        load(sprintf('%s/x_ob_unseen.%s.mat', save_dir, save_id), 'x_unseen', 'tempsaveseed');
        if exist('tempsaveseed', 'var'),
           rng(tempsaveseed);
           clear tempsaveseed
        end
        % Note that if x_ob_unseen was saved without saving the seeds at that point,
        % then the run will still happen, and still be reproducible, but not be the
        % same as if it had been completely redone from scratch.
    else
       trained_smpl_file = sample_Stu_mix(training_data, primary_prior_file, primary_smpl_num, ...
                                          save_dir, save_id, primary_ctrl_file, [], analysis_seed);

       x_unseen = x_ob(:, is_unseen);
       tempsaveseed = rng;

       % Save unseen x for testing.
       save(sprintf('%s/x_ob_unseen.%s.mat', save_dir, save_id), 'x_unseen', 'tempsaveseed', '-v6');
    end

    tmp_j = calc_j(trained_smpl_file, x_unseen, lambda, Nskip); % Again, data not passed to file.
    j_vec((it_N - 1) * N_unseen + 1: it_N * N_unseen) = tmp_j;

    fprintf('Saving %s for it_N = %d of %d...\n', ...
            sprintf('%s/tempsave.%s.mat', save_dir, primary_id), ...
            it_N, Nsplit);
    tempsaveseed = rng;
    save(sprintf('%s/tempsave.%s.mat', save_dir, primary_id));
    fprintf('...tempsave finished.\n');
    clear tempsaveseed

end

if ~isempty(analysis_seed) % Restore previous settings if rng was set using given analysis seed.
     rng(previousSettings);
end

analysis_seed = tmp_analysis_seed; % Saved as the seed used in this run.

x_ob = j_vec;
is_synth = 0;
ASI_data_file = sprintf('%s/ASI_data.%s.mat', save_dir, ASI_id); 
save(ASI_data_file, '-v6');

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
