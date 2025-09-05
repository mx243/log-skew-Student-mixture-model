function plot_curves(all_curves, plot_num, plotx, curve025, curve975, curve_mean, ...
                     is_synth, true_curve, cp1, cp2, titl)
% Usage: plot_curves(all_curves, plot_num, plotx, curve025, curve975, curve_mean, ...
%                    is_synth, true_curve, cp1, cp2, titl)
% This function makes two plots. The first one contains all curves
% specified by plotx and the rows of all_curves, and optionally true_curve.
% The second one contains curve025, curve975, curve_mean and optionally
% true_curve.
% (Will add Kaplan_Meier)

% Change Log:
%
%     1.1          03:oct:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_curves.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

figure; % All samples.
for num = 1 : plot_num
    hold on;
    plot(plotx, all_curves(num, :), cp1);
end
if is_synth
    hold on;
    plot(plotx, true_curve, 'g-', 'LineWidth', 1.5);
end
title([titl, ': all samples']);

figure; % Mean and centiles.
if is_synth
    hold on;
    plot(plotx, true_curve, 'g-', 'LineWidth', 1.5);
end
hold on;
plot(plotx, curve025, cp2);
hold on;
plot(plotx, curve975, cp2);
hold on;
plot(plotx, curve_mean, cp1);
title([titl, ': mean and centiles']);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
