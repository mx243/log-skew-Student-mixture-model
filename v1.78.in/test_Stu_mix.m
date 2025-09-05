function test_Stu_mix(data_file, prior_file, save_dir, save_id, ctrl_file, Enables, analysis_seed)
% Usage: test_Stu_mix(data_file, prior_file, save_dir, save_id, ctrl_file, Enables, analysis_seed)
% This function reads synthetic data from data_file, calls 'sample_Stu_mix'
% with prior_file and ctrl_file to get [Samples, Inter_smpls, I_thermo] and
% save these results at save_dir with save_id.
% It can reproduce a set of results when given analysis_seed.
% The function presents the results (p, m, \mu, S, \nu) in two different
% ways: e.g. for p, 
% 1. Plot the true values of p_1, ..., p_min(7, nonempty_c_num) as
%    horizontal lines, then for each sample, plot p_1, ..., 
%    p_min(7, nonempty_c_num) as points.
% 2. Find k_1, ..., k_min(7, nonempty_c_num) with different colours. Then
%    plot both the true values and the sample values of p_k_1, ..., 
%    p_k_min(7, nonempty_c_num) as lines.
% The values of axes above \phi in the model diagram are presented as
% usual. For alpha_1, alpha_2, ..., only plot the first 7.
% For axes like mu_mu that have dimension > 1, plot the first 5 entries in
% separate figures.

% Change Log:
%
%     1.1          04:sep:24    mx243      First version.
%     1.20         09:sep:24    rfs34      Warning added that this function is now redundant.
%     1.26         13:sep:24    rfs34      Despite that updated to take account of typ_m.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.test_Stu_mix.m 1.78 25/06/09 17:25:32 ';

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

N = size(Inter_smpls.nonempty_c_num, 2);
K = size(Inter_smpls.alpha, 1);
D = Priors.D;
smpl_cnt_vec = 1 : N;

cols = 'rgbcmyk';
Ncols = length(cols);
Ndim = 5;

% Plot alpha.
if Enables.alpha
    clf;
    for k = 1 : min(K, Ncols)
        hold on;
        plot(smpl_cnt_vec, Inter_smpls.alpha(k, :), [cols(k), '-']);
        % Plot true value in the same colour but thicker line.
        plot([1, N], Args_true.alpha(k) * ones(1, 2), [cols(k), '-'], 'LineWidth', 1.5);
        hold off;
    end
    title('\alpha');
    pause;
end

% Plot R_S_mu.
if Enables.R_S_mu
    clf;
    nplot = 0;
    for row = 1 : D
        for col = 1 : D
            nplot = nplot + 1;
            if nplot > Ncols
                break;
            end
            hold on;
            plot(smpl_cnt_vec, squeeze(Inter_smpls.R_S_mu(row, col, :)), [cols(nplot), '-']);
            plot([1, N], Args_true.R_S_mu(row, col) * ones(1, 2), [cols(nplot), '-'], 'LineWidth', 1.5);
            hold off;
        end
    end
    title('R_{S_{\mu}}');
    pause;
end

% Plot S_mu.
if Enables.S_mu
    clf;
    nplot = 0;
    for row = 1 : D
        for col = 1 : D
            nplot = nplot + 1;
            if nplot > Ncols
                break;
            end
            hold on;
            plot(smpl_cnt_vec, squeeze(Inter_smpls.S_mu(row, col, :)), [cols(nplot), '-']);
            plot([1, N], Args_true.S_mu(row, col) * ones(1, 2), [cols(nplot), '-'], 'LineWidth', 1.5);
            hold off;
        end
    end
    title('S_{\mu}');
    pause;
end

% Plot mu_mu.
if Enables.mu_mu
    clf;
    for k = 1 : min(D, Ncols)
        hold on;
        plot(smpl_cnt_vec, Inter_smpls.mu_mu(k, :), [cols(k), '-']);
        plot([1, N], Args_true.mu_mu(k) * ones(1, 2), [cols(k), '-'], 'LineWidth', 1.5);
        hold off;
    end
    title('\mu_{\mu}');
    pause;
end

% Plot m_S.
if Enables.m_S
    clf;
    hold on;
    plot(smpl_cnt_vec, Inter_smpls.m_S(:), 'b-');
    plot([1, N], Args_true.m_S * ones(1, 2), 'g-', 'LineWidth', 1.5);
    hold off;
    title('m_S');
    pause;
end

% Plot R_S.
if Enables.R_S
    clf;
    nplot = 0;
    for row = 1 : D
        for col = 1 : D
            nplot = nplot + 1;
            if nplot > Ncols
                break;
            end
            hold on;
            plot(smpl_cnt_vec, squeeze(Inter_smpls.R_S(row, col, :)), [cols(nplot), '-']);
            plot([1, N], Args_true.R_S(row, col) * ones(1, 2), [cols(nplot), '-'], 'LineWidth', 1.5);
            hold off;
        end
    end
    title('R_S');
    pause;
end

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
