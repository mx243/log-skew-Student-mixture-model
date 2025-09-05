function plot_smpls(Samples, LN, N, is_prior, Args_true, is_synth, x1_bar, T_max, Npt)
% Usage: plot_smpls(Samples, LN, N, is_prior, Args_true, is_synth, x1_bar, T_max, Npt)
% This function plots the survival probability, lifetime and hazard rate at
% one x1_bar given by each sample of distribution arguments, and their 
% respective mean and centiles over all samples. If the samples are
% posterior samples obtained from synthetic data, the true distributions
% will also be plotted.

% Change Log:
%
%     1.1          01:oct:24    mx243      First version.
%     1.47         04:oct:24    rfs34      Layout only changed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_smpls.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

D = size(x1_bar, 1) + 1;
Phi = Samples.phi(LN : N);
C = Samples.C(LN : N);
smpl_num = N - LN + 1; % Number of samples we extracted at t_SA = 1.

p_dis = 0.15; % Discarding the first 15% of samples.
start = ceil(p_dis * smpl_num);
plot_num = smpl_num - start + 1;

if is_prior % Color and pattern for the samples.
    cp1 = 'm-';
    cp2 = 'm--';
    is_synth = 0; % Whether true arguments exist is irrelevant when the samples came from the priors.
else
    cp1 = 'b-';
    cp2 = 'b--';
end

% x and y for the samples.
plotx = 0 : T_max / Npt : T_max;
log_lifetime = -Inf(plot_num, Npt + 1);
survival = NaN(plot_num, Npt + 1);
hazardrt = NaN(plot_num, Npt + 1);
tmp_log_life = NaN(1, Npt + 1);

% tmpx will be used as an argument of log_density_skewStu.
tmpx = NaN(D, Npt + 1);
tmpx(1, :) = log(plotx);
tmpx(2 : D, :) = x1_bar(:, ones(1, Npt + 1));

T = eye(D);
T = T(2 : D, :);
cnt = 0;
for num = start : smpl_num 
    cnt = cnt + 1;
    for it_c = 1 : C(num)
        mu = Phi(num).mu(:, it_c);
        S = Phi(num).S(:, :, it_c);
        m = Phi(num).m(it_c);
        nu = Phi(num).nu(:, it_c);
        mu_bar = mu(2 : D);
        S_bar = (T * (S \ T')) ^ -1;
        nu_bar = nu(2 : D);

        % Pdf of x1_bar evaluated at the given x1_bar. Note r = m in the model for primary data analysis.
        norm_const = log_density_skewStu(x1_bar, mu_bar, S_bar, m, m, nu_bar, D - 1); 

        tmp_log_life = log_density_skewStu(tmpx, mu, S, m, m, nu, D) - norm_const; % Normalise to conditional prob.
        tmp_log_life = tmp_log_life + log(Phi(num).p(it_c)); % Weighted average wrt p.

        log_lifetime(cnt, :) = -Eadd(-log_lifetime(cnt, : ), -tmp_log_life); % Here we have the log of the lifetime pdf.
    end
    log_lifetime(cnt, :) = log_lifetime(cnt, :) - log(plotx); % Factor from the chain rule.
    log_lifetime(cnt, 1) = -Inf; % Otherwise would be NaN.
end

lifetime = exp(log_lifetime);

% Use numerical integration which is much faster than doing Monte-Carlo integration over alpha.
survival = 1 - cumtrapz(plotx, lifetime, 2); 

% I don't know why plotx has to be the first entry, while it's used to
%  calculate numerical gradient along the second entry of survival.
hazardrt = - gradient(survival, plotx, 1) ./ survival;

LB = ceil(0.025 * plot_num);
UB = floor(0.975 * plot_num);

% Mean and centiles of lifetime.
tmp_lifetime = sort(lifetime, 1);
lifetime025 = tmp_lifetime(LB, :);
lifetime975 = tmp_lifetime(UB, :);
lifetime_mean = exp(- Esum(-log_lifetime, 1) - log(plot_num));

% Mean and centiles of survival.
tmp_survival = sort(survival, 1);
survival025 = tmp_survival(LB, :);
survival975 = tmp_survival(UB, :);
survival_mean = mean(survival, 1);

% Mean and centiles of hazardrt.
tmp_hazardrt = sort(hazardrt, 1);
hazardrt025 = tmp_hazardrt(LB, :);
hazardrt975 = tmp_hazardrt(UB, :);
hazardrt_mean = mean(hazardrt, 1);

lifetime_true = NaN(1, Npt + 1); % These will still be fed into plot_curves if is_synth = 0, for easier coding.
survival_true = NaN(1, Npt + 1);
hazardrt_true = NaN(1, Npt + 1);
if is_synth % Compute the true distributions if the data is synthetic.
    log_lifetime_true = -Inf(1, Npt + 1);
    for it_c = 1 : Args_true.phi.C
        mu = Args_true.phi.mu(:, it_c);
        S = Args_true.phi.S(:, :, it_c);
        m = Args_true.phi.m(it_c);
        nu = Args_true.phi.nu(:, it_c);
        mu_bar = mu(2 : D);
        S_bar = (T * (S \ T')) ^ -1;
        nu_bar = nu(2 : D);

        % Pdf of x1_bar evaluated at the given x1_bar. Note r = m in the model for primary data analysis.
        norm_const = log_density_skewStu(x1_bar, mu_bar, S_bar, m, m, nu_bar, D - 1); 

        tmp_log_life = log_density_skewStu(tmpx, mu, S, m, m, nu, D) - norm_const; % Normalise to conditional prob.
        tmp_log_life = tmp_log_life + log(Args_true.phi.p(it_c)); % Weighted average wrt p.

        log_lifetime_true = -Eadd(-log_lifetime_true, -tmp_log_life); % Here we have the log of the lifetime pdf.
    end
    log_lifetime_true = log_lifetime_true - log(plotx); % Factor from the chain rule.
    log_lifetime_true(1) = -Inf; % Otherwise would be NaN.

    lifetime_true = exp(log_lifetime_true);
    survival_true = 1 - cumtrapz(plotx, lifetime_true, 2);
    hazardrt_true = - gradient(survival_true, plotx, 1) ./ survival_true;
end

close all;

% Plot lifetime.
plot_curves(lifetime, plot_num, plotx, lifetime025, lifetime975, lifetime_mean, ...
            is_synth, lifetime_true, cp1, cp2, 'lifetime distribution');
% Plot survival.
plot_curves(survival, plot_num, plotx, survival025, survival975, survival_mean, ...
            is_synth, survival_true, cp1, cp2, 'survival probability');
% Plot hazardrt.
plot_curves(hazardrt, plot_num, plotx, hazardrt025, hazardrt975, hazardrt_mean, ...
            is_synth, hazardrt_true, cp1, cp2, 'hazard rate');

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
