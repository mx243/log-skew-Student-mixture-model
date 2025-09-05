function plot_post_smpls(Samples, Args_true, x1_bar, x_ob, T_max, Npt, mode, if_centiles, typ, usealpha, Nskip)
% Usage: plot_post_smpls(Samples, Args_true, x1_bar, x_ob, T_max, Npt, mode, if_centiles, typ, usealpha, Nskip)
% This function plots one or two sets of posterior samples. When plotting one,
% the curves are in blue. The true values (if exist) are plotted in green.
% If there's a second set, the first is plotted in red and the second in green. 
% We do not envisage plotting two sets of samples and the true ones on the same plot.
% if_centiles = 0: plot the samples themselves
% if_centiles = 1: plot the mean and centiles
% if_centiles = 2: plot both
% When there're two sets of samples, they must be obtained from the same
% set of training data, i.e. Args_true (if exist) and x_ob (hence the K-M plot) 
% must be the same.
% The function can only operate on one set of x1_bar's and one mode.
% Samples should be a cell array of size {1, Ngroups}.
% typ is the proGamma type used for the prior on m, i.e. Priors.typ_m.
% This is normally 2 for the primary data analysis, but is 1 if distribution of ASI is being modelled.
%
% usealpha = 0: survival data is calculated by doing numerical integration
%               to lifetime pdf.
% usealpha > 0: survival data is calculated by doing Monte-Carlo integration over alpha using usealpha
%               samples of alpha.
% usealpha < 0: survival data is calculated by doing numerical integration over alpha using abs(usealpha)
%               uniformly log-spaced values of alpha.
% Nskip allows only every Nskipth sample to be used.

% Change Log:
%
%     1.1          11:oct:24    mx243      First version.
%     1.52         14:oct:24    rfs34      Colour conventions changed to match previous ones;
%                                          typ required as a parameter; comments corrected.
%     1.54         18:oct:24    mx243      Added usealpha, and suppress the
%                                          K-M plot when x_ob is empty.
%     1.57         20:oct:24    mx243      Completed the usealpha option
%                                          for the true values.
%     1.58         23:oct:24    mx243      Plot the truth in the same mode as the samples.
%                                          Added extra plots illustrating the comparison of 
%                                          survival probability and hazard rate
%                                          between two sets of samples.
%     1.60         26:oct:24    rfs34      Fixed bug when Args_true present but actually empty;
%                                          added Nskip to parameters.
%     1.71         25:jan:25    rfs34      Comments only changed.
%     1.75         19:mar:25    mx243      Adapt to changes in get_probs, etc.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_post_smpls.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 11 || isempty(Nskip),
   Nskip = 1;
end

if mode == 'A'
    subtitl = 'patients'; % Each curve corresp to a different patient.
else
    subtitl = 'samples'; % Each curve corresp to a different sample.
end

if ~isempty(x_ob)
    x_ob1 = x_ob(1, :);
    is_censored = x_ob1 < 0;
    x_ob1 = abs(x_ob1);
end

sz = size(Samples, 2);
if sz == 2,
    % Changed to match Tommy's color convention.
    color_all = {'r-', 'g-'};
    color_cent = {'r--', 'g--'};

    % Indicate that red is group 1, green is group 2.
    subtitl = [subtitl, ', R1G2'];

    if ~isempty(Args_true),
       error('Not expecting two sets of plots and also the true values');
    end
elseif sz == 1,
    color_all = {'b-'};
    color_cent = {'b--'};
end

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

    % Distribution of x, ddf P(x1 > d, x1bar) and marginal distribution of x1bar:
    [log_Px, log_ddf, log_Px1bar] = get_probs(Phi, C, plotx, typ, usealpha, x1_bar); 

    if ~usealpha
        [lifetime{cnt}, survival{cnt}, hazardrt{cnt}, MClifetime{cnt}, MCsurvival{cnt}, MChazardrt{cnt}] = ...
        get_prediction(log_Px, log_Px1bar, plotx, mode, if_centiles, usealpha);
    else
        [lifetime{cnt}, survival{cnt}, hazardrt{cnt}, MClifetime{cnt}, MCsurvival{cnt}, MChazardrt{cnt}] = ...
        get_prediction_use_alpha(log_Px, log_ddf, log_Px1bar, plotx, mode, if_centiles, usealpha);
    end

end

if ~isempty(Args_true) && Args_true.phi.C > 0,

    % Joint and marginal distributions:
    [log_Px, log_ddf, log_Px1bar] = get_probs(Args_true.phi, Args_true.phi.C, plotx, typ, usealpha, x1_bar); 

    % Use mode A, i.e. the prediction for each patient obtained by
    % 'averaging' (see get_prediction) over all samples of theta. 
    % In this case, there's only one theta which is the true value.
    %
    % Now switched to the same mode as for the samples, for easier
    % comparison.
    if ~usealpha
        [lft_true, svl_true, hzd_true, MClft_true, MCsvl_true, MChzd_true] = ...
        get_prediction(log_Px, log_Px1bar, plotx, mode, if_centiles, usealpha);
    else
        [lft_true, svl_true, hzd_true, MClft_true, MCsvl_true, MChzd_true] = ...
        get_prediction_use_alpha(log_Px, log_ddf, log_Px1bar, plotx, mode, if_centiles, usealpha);
    end

    sz = sz + 1; % One more set of values to plot.

    color_all{sz} = 'g-';
    color_cent{sz} = 'g--';

    % Expand to include true values.
    lifetime{sz} = lft_true;
    survival{sz} = svl_true;
    hazardrt{sz} = hzd_true;
    MClifetime{sz} = MClft_true;
    MCsurvival{sz} = MCsvl_true;
    MChazardrt{sz} = MChzd_true;
end

% close all;

if if_centiles == 0 || if_centiles == 2

    ord = [];
    for cnt = 1 : sz
        ord = [ord, cnt * ones(1, size(lifetime{cnt}, 1))];
    end
    ord = ord(randperm(size(ord, 2))); % Randomise order of plottings.

    titl = ['Post samples of different ', subtitl];

    plot_all_curves(plotx, lifetime, sz, ord, color_all, titl, 'Time(days)', 'Lifetime pdf');

    plot_all_curves(plotx, survival, sz, ord, color_all, titl, 'Time(days)', 'Survival probability');
    hold on;
    if ~isempty(x_ob)
        plot_survival_data(x_ob1', is_censored', T_max, 'k-');
    end

    plot_all_curves(plotx, hazardrt, sz, ord, color_all, titl, 'Time(days)', 'Hazard rate');

end

if if_centiles == 1 || if_centiles == 2

    if mode == 'C'
        subtitl = 'patients and samples';
        if sz == 2
            subtitl = [subtitl, ', R1G2'];
        end
    end

    titl = ['Mean and centiles of Post samples of different ', subtitl];

    plot_mean_and_centiles(plotx, MClifetime, sz, color_all, color_cent, titl, 'Time(days)', 'Lifetime pdf');

    plot_mean_and_centiles(plotx, MCsurvival, sz, color_all, color_cent, titl, 'Time(days)', 'Survival probability');
    hold on;
    if ~isempty(x_ob)
        plot_survival_data(x_ob1', is_censored', T_max, 'k-');
    end
    
    plot_mean_and_centiles(plotx, MChazardrt, sz, color_all, color_cent, titl, 'Time(days)', 'Hazard rate');

end

if sz == 2 % (We are also interested in where the truth is in the posterior.)

    plot_comparison_prob(plotx, survival, sz, 'Survival probability');

    plot_comparison_prob(plotx, hazardrt, sz, 'Hazard rate');

end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
