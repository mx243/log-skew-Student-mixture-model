function plot_prior_smpls(Samples, T_max, Npt, if_centiles, typ, usealpha, Nskip)
% Usage: plot_prior_smpls(Samples, T_max, Npt, if_centiles, typ, usealpha, Nskip)
% This function plots one or two sets of prior samples. When plotting one,
% the curves are in magenta. If there's a second set, they're plotted in
% cyan. 
% if_centiles = 0: plot the samples themselves
% if_centiles = 1: plot the mean and centiles
% if_centiles = 2: plot both
% typ is the proGamma type used for the prior on m, i.e. Priors.typ_m.
% This is normally 2 for the primary data analysis, but is 1 if distribution of ASI is being modelled.
%
% usealpha = 0: survival data is calculated by doing numerical integration
%               to lifetime pdf.
% usealpha > 0: survival data is calculated by doing Monte-Carlo integration over alpha,
%               using usealpha samples of alpha.
% usealpha < 0: survival data is calculated by doing numerical integration over alpha,
%               using abs(usealpha) uniformly log-spaced values of alpha.
% Nskip allows only every Nskipth sample to be used.

% Change Log:
%
%     1.1          08:oct:24    mx243      First version.
%     1.52         14:oct:24    rfs34      typ required as a parameter.
%     1.54         17:oct:24    mx243      Added usealpha.
%     1.60         26:oct:24    rfs34      Added Nskip option.
%     1.71         25:jan:25    rfs34      Comments only changed.
%     1.75         19:mar:25    mx243      Adapt to changes in get_probs, etc.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_prior_smpls.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 7 || isempty(Nskip),
   Nskip = 1;
end

sz = size(Samples, 2);
color_all = {'m-', 'c-'};
color_cent = {'m--', 'c--'};
plotx = 0 : T_max / Npt : T_max;

lifetime = cell(1, sz);
survival = cell(1, sz);
hazardrt = cell(1, sz);
MClifetime = cell(1, sz);
MCsurvival = cell(1, sz);
MChazardrt = cell(1, sz);

for cnt = 1 : sz

    LN0 = Samples{cnt}.LN0;
    LN1 = Samples{cnt}.LN1;

    p_dis = 0.15; % Discarding the first 15% of samples.
    start = LN0 + ceil(p_dis * (LN1 - LN0 + 1));

    Phi = Samples{cnt}.phi(start : Nskip : LN1);
    C = Samples{cnt}.C(start : Nskip : LN1);

    [log_Px, log_ddf, log_Px1bar] = get_probs(Phi, C, plotx, typ, usealpha); % Marginal distr of x1.

    if ~usealpha
        % Use type B here, i.e. for each sample, average across all patients. But
        % here we have no x1_bar, where we put Npatient = 1, hence no averaging
        % going on.
        [lifetime{cnt}, survival{cnt}, hazardrt{cnt}, MClifetime{cnt}, MCsurvival{cnt}, MChazardrt{cnt}] = ...
        get_prediction(log_Px, log_Px1bar, plotx, 'B', if_centiles, usealpha);
    else
        [lifetime{cnt}, survival{cnt}, hazardrt{cnt}, MClifetime{cnt}, MCsurvival{cnt}, MChazardrt{cnt}] = ...
        get_prediction_use_alpha(log_Px, log_ddf, log_Px1bar, plotx, 'B', if_centiles, usealpha);
    end
end

% close all;

if if_centiles == 0 || if_centiles == 2

    if sz == 2
        ord = [ones(1, size(lifetime{1}, 1)), 2 * ones(1, size(lifetime{2}, 1))];
        ord = ord(randperm(size(ord, 2))); % Randomise order of plottings.
    else
        ord = ones(1, size(lifetime{1}, 1));
    end

    plot_all_curves(plotx, lifetime, sz, ord, color_all, 'Prior samples', 'Time(days)', 'Lifetime pdf');

    plot_all_curves(plotx, survival, sz, ord, color_all, 'Prior samples', 'Time(days)', 'Survival probability');
    ylim([-0.01, 1.01]);

    plot_all_curves(plotx, hazardrt, sz, ord, color_all, 'Prior samples', 'Time(days)', 'Hazard rate');
    ylim([0, 0.3]);

end

if if_centiles == 1 || if_centiles == 2
    
    plot_mean_and_centiles(plotx, MClifetime, sz, color_all, color_cent, 'Prior mean and centiles', 'Time(days)', 'Lifetime pdf');

    plot_mean_and_centiles(plotx, MCsurvival, sz, color_all, color_cent, 'Prior mean and centiles', 'Time(days)', 'Survival probability');
    ylim([-0.01, 1.01]);

    plot_mean_and_centiles(plotx, MChazardrt, sz, color_all, color_cent, 'Prior mean and centiles', 'Time(days)', 'Hazard rate');
    ylim([0, 0.9]);

end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
