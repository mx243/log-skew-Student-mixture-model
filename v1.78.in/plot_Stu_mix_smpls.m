function plot_Stu_mix_smpls(smpl_file)
% Usage: plot_Stu_mix_smpls(smpl_file)
% This function plots the samples from smpl_file.
% It presents the results (p, m, \mu, S, \nu) in two different
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
%     1.1          05:sep:24    mx243      First version.
%     1.20         09:sep:24    rfs34      Includes all variables including x1; Latex titles
%                                          etc now work; each plot in a new figure so that figures
%                                          can easily be compared; non-empty cluster number suppressed. 
%     1.54         17:oct:24    rfs34      Now fails gracefully if Args_true not fully set up.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_Stu_mix_smpls.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

load(smpl_file);

N = size(Samples.nonempty_c_num, 2);
K = size(Samples.alpha, 1);
D = Priors.D;
smpl_cnt_vec = 1 : N;

cols = 'rgbcmyk';
Ncols = length(cols);
Ndim = 5;

% Plot alpha.
if Enables.alpha
    figure;
    clf;
    for k = 1 : min(K, Ncols)
        hold on;
        plot(smpl_cnt_vec, Samples.alpha(k, :), [cols(k), '-']);
        % Plot true value in the same colour but thicker line.
        plot([1, N], Args_true.alpha(k) * ones(1, 2), [cols(k), '-'], 'LineWidth', 1.5);
        hold off;
    end
    title('$\alpha$', 'Interpreter', 'Latex');
    keyboard;
end

if Enables.x1,
   figure;
   clf;
   for k = 1 : min(K, Ncols),
      hold on;
      plot(smpl_cnt_vec, Samples.x1(k, :), [cols(k), '-']);
      % Plot true value in same colour but thicker line.
      plot([1, N], Args_true.x1(ones(2, 1), k).', [cols(k), '-x'], 'LineWidth', 1.5, 'MarkerSize', 10);
      hold off;
   end
   title('x1');
   keyboard;
end

% Plot R_S_mu.
if Enables.R_S_mu
    figure;
    clf;
    nplot = 0;
    for row = 1 : D
        for col = 1 : D
            nplot = nplot + 1;
            if nplot > Ncols
                break;
            end
            hold on;
            plot(smpl_cnt_vec, squeeze(Samples.R_S_mu(row, col, :)), [cols(nplot), '-']);
            plot([1, N], Args_true.R_S_mu(row, col) * ones(1, 2), [cols(nplot), '-'], 'LineWidth', 1.5);
            hold off;
        end
    end
    title('$R_{S_{\mu}}$', 'Interpreter', 'Latex');
    keyboard;
end

% Plot S_mu.
if Enables.S_mu
    figure;
    clf;
    nplot = 0;
    for row = 1 : D
        for col = 1 : D
            nplot = nplot + 1;
            if nplot > Ncols
                break;
            end
            hold on;
            plot(smpl_cnt_vec, squeeze(Samples.S_mu(row, col, :)), [cols(nplot), '-']);
            plot([1, N], Args_true.S_mu(row, col) * ones(1, 2), [cols(nplot), '-'], 'LineWidth', 1.5);
            hold off;
        end
    end
    title('$S_{\mu}$', 'Interpreter', 'Latex');
    keyboard;
end

% Plot mu_mu.
if Enables.mu_mu
    figure;
    clf;
    for k = 1 : min(D, Ncols)
        hold on;
        plot(smpl_cnt_vec, Samples.mu_mu(k, :), [cols(k), '-']);
        plot([1, N], Args_true.mu_mu(k) * ones(1, 2), [cols(k), '-'], 'LineWidth', 1.5);
        hold off;
    end
    title('$\mu_{\mu}$', 'Interpreter', 'Latex');
    keyboard;
end

% Plot m_S.
if Enables.m_S
    figure;
    clf;
    hold on;
    plot(smpl_cnt_vec, Samples.m_S(:), 'b-');
    plot([1, N], Args_true.m_S * ones(1, 2), 'g-', 'LineWidth', 1.5);
    hold off;
    title('m_S');
    keyboard;
end

% Plot R_S.
if Enables.R_S
    figure;
    clf;
    nplot = 0;
    for row = 1 : D
        for col = 1 : D
            nplot = nplot + 1;
            if nplot > Ncols
                break;
            end
            hold on;
            plot(smpl_cnt_vec, squeeze(Samples.R_S(row, col, :)), [cols(nplot), '-']);
            plot([1, N], Args_true.R_S(row, col) * ones(1, 2), [cols(nplot), '-'], 'LineWidth', 1.5);
            hold off;
        end
    end
    title('R_S');
    keyboard;
end

if Enables.phi.C
    figure;
    clf;
    hold on;
    plot(smpl_cnt_vec, Samples.C(:), 'b-', 'LineWidth', 2);
    plot([1, N], Args_true.phi.C * ones(1, 2), 'g-x', 'LineWidth', 2, 'MarkerSize', 10);
    % plot(smpl_cnt_vec, Samples.nonempty_c_num(:), 'r-x', 'LineWidth', 1);
    % plot([1, N], Args_true.nonempty_c_num * ones(1, 2), 'm', 'LineWidth', 1);
    hold off;
    % title('C and nonempty c num');
    title('C');
    keyboard;
end

if Args_true.phi.C == 0,
   fprintf('The remaining code in this function, for plotting samples of phi, will not work without Args_true being set up properly.\n');
   fprintf('Try using plotphi.m instead.\n');
   return;
end

% Find some k's with distinct c(k).
plot_k_num = min(Ncols, max([0; Args_true.c(:)]));
plot_k = NaN(1, plot_k_num);
seen_c = zeros(1, max([0; Args_true.c(:)]));
k_cnt = 1;
for k = 1 : K
    if k_cnt > plot_k_num
        break;
    end
    if ~seen_c(Args_true.c(k))
        seen_c(Args_true.c(k)) = 1;
        plot_k(k_cnt) = k;
        k_cnt = k_cnt + 1;
    end
end

if Enables.phi.p
    %{
    figure;
    clf;
    hold on;
    for t = 1 : 1000
        for it_c = 1 : min(Inter_smpls.phi(floor(t * N / 1000)).C, Ncols)
            % hold on;
            plot(floor(t * N / 1000), Inter_smpls.phi(floor(t * N / 1000)).p(it_c), [cols(it_c), '.'], 'MarkerSize', 10);
        end
    end
    for it_c = 1 : min(Args_true.phi.C, Ncols)
        % hold on;
        plot([1, N], Args_true.phi.p(it_c) * ones(1, 2), [cols(it_c), '-'], 'LineWidth', 1.5);
    end
    hold off;
    title('p: plot 1');

    pause;
    %}
    figure;
    clf;
    p_c_k = NaN(plot_k_num, N);
    hold on;
    for it_pos = 1 : plot_k_num
        for it_N = 1 : N
            % if(Samples.c(plot_k(it_pos)) > 6), keyboard; end;
            p_c_k(it_pos, it_N) = Samples.phi(it_N).p(Samples.c(plot_k(it_pos), it_N));
        end
        plot(smpl_cnt_vec, p_c_k(it_pos, :), [cols(it_pos), '-']);
        plot([1, N], Args_true.phi.p(Args_true.c(plot_k(it_pos))) * ones(1, 2), ...
             [cols(it_pos), '-x'], 'LineWidth', 1.5, 'MarkerSize', 10);
    end
    hold off;
    title('p: plot 2');
    keyboard;
end

if Enables.phi.m
    %{
    figure;
    clf;
    hold on;
    for t = 1 : 1000
        for it_c = 1 : min(Inter_smpls.phi(floor(t * N / 1000)).C, Ncols)
            % hold on;
            plot(floor(t * N / 1000), Inter_smpls.phi(floor(t * N / 1000)).m(it_c), [cols(it_c), '.'], 'MarkerSize', 10);
        end
    end
    for it_c = 1 : min(Args_true.phi.C, Ncols)
        % hold on;
        plot([1, N], Args_true.phi.m(it_c) * ones(1, 2), [cols(it_c), '-'], 'LineWidth', 1.5);
    end
    hold off;
    title('m: plot 1');

    pause;
%}
    figure;
    clf;
    m_c_k = NaN(plot_k_num, N);
    hold on;
    for it_pos = 1 : plot_k_num
        % if it_pos ~= 3, continue; end
        for it_N = 1 : N
            m_c_k(it_pos, it_N) = Samples.phi(it_N).m(Samples.c(plot_k(it_pos), it_N)); 
        end
        plot(smpl_cnt_vec, m_c_k(it_pos, :), [cols(it_pos), '-']);
        plot([1, N], Args_true.phi.m(Args_true.c(plot_k(it_pos))) * ones(1, 2), [cols(it_pos), '-'], 'LineWidth', 1.5);
    end
    hold off;
    title('m: plot 2');
    keyboard;
end

if Enables.phi.S
    ax_cnt = 0;
    for row = 1 : D
        for col = 1 : D
            ax_cnt = ax_cnt + 1;
            if ax_cnt > 5
                break;
            end
            figure;
            clf;
            S_c_k = NaN(plot_k_num, N);
            hold on;
            for it_pos = 1 : plot_k_num
                % if it_pos ~= 3, continue; end
                for it_N = 1 : N
                    S_c_k(it_pos, it_N) = Samples.phi(it_N).S(row, col, Samples.c(plot_k(it_pos), it_N)); 
                end
                plot(smpl_cnt_vec, S_c_k(it_pos, :), [cols(it_pos), '-']);
                plot([1, N], Args_true.phi.S(row, col, Args_true.c(plot_k(it_pos))) * ones(1, 2), [cols(it_pos), '-'], 'LineWidth', 1.5);
            end
            hold off;
            title(sprintf('S(%d, %d): plot 2', row, col));
            keyboard;
        end
    end
end

    %{
    figure;
    clf;
    hold on;
    for t = 1 : 1000
        for it_c = 1 : min(Inter_smpls.phi(floor(t * N / 1000)).C, Ncols)
            % hold on;
            plot(floor(t * N / 1000), Inter_smpls.phi(floor(t * N / 1000)).m(it_c), [cols(it_c), '.'], 'MarkerSize', 10);
        end
    end
    for it_c = 1 : min(Args_true.phi.C, Ncols)
        % hold on;
        plot([1, N], Args_true.phi.m(it_c) * ones(1, 2), [cols(it_c), '-'], 'LineWidth', 1.5);
    end
    hold off;
    title('m: plot 1');

    pause;
%}

if Enables.phi.mu
    ax_cnt = 0;
    for row = 1 : D
        ax_cnt = ax_cnt + 1;
        if ax_cnt > 5
            break;
        end
        figure;
        clf;
        mu_c_k = NaN(plot_k_num, N);
        hold on;
        for it_pos = 1 : plot_k_num
            % if it_pos ~= 3, continue; end
            for it_N = 1 : N
                mu_c_k(it_pos, it_N) = Samples.phi(it_N).mu(row, Samples.c(plot_k(it_pos), it_N)); 
            end
            plot(smpl_cnt_vec, mu_c_k(it_pos, :), [cols(it_pos), '-']);
            plot([1, N], Args_true.phi.mu(row, Args_true.c(plot_k(it_pos))) * ones(1, 2), [cols(it_pos), '-'], 'LineWidth', 1.5);
        end
        hold off;
        title(sprintf('mu(%d): plot 2', row));
        keyboard;
    end
end

if Enables.phi.nu
    ax_cnt = 0;
    for row = 1 : D
        ax_cnt = ax_cnt + 1;
        if ax_cnt > 5
            break;
        end
        figure;
        clf;
        nu_c_k = NaN(plot_k_num, N);
        hold on;
        for it_pos = 1 : plot_k_num
            % if it_pos ~= 3, continue; end
            for it_N = 1 : N
                nu_c_k(it_pos, it_N) = Samples.phi(it_N).nu(row, Samples.c(plot_k(it_pos), it_N)); 
            end
            plot(smpl_cnt_vec, nu_c_k(it_pos, :), [cols(it_pos), '-']);
            plot([1, N], Args_true.phi.nu(row, Args_true.c(plot_k(it_pos))) * ones(1, 2), [cols(it_pos), '-'], 'LineWidth', 1.5);
        end
        hold off;
        title(sprintf('$\\nu_{%d}$: plot 2', row), 'Interpreter', 'Latex');
        keyboard;
    end
end

return;

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
