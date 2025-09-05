function plot_mean_and_centiles(plotx, MC, Ngroups, col1, col2, titl, xlab, ylab)
% Usage: plot_mean_and_centiles(plotx, MC, Ngroups, col1, col2, titl, xlab, ylab)
% This function plots the (possibly two sets of) mean and centiles of
% samples.
% plotx should be a vector of size [1, Nx].
% MC should be a cell array of structures of size {1, Ngroups},
%  with MC{ngroup}.mean being a vector of size [1, Nx] giving the mean
%  and similarly for MC{ngroup}.p025 and MC{ngroup}.p975 giving the centiles.
% col1 should be the colour and linestyle string for plotting the mean,
%  and col2 that for plotting the centiles.
% titl, xlab, and ylab are respectively the desired title, xlabel, and ylabel for the plot.

% Change Log:
%
%     1.1          10:oct:24    mx243      First version.
%     1.52         14:oct:24    rfs34      Comments extended and sz renamed to Ngroups, similarly cnt.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_mean_and_centiles.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

figure;
clf;
for ngroup = 1 : Ngroups
    hold on;
    plot(plotx, MC{ngroup}.mean, col1{ngroup}, 'LineWidth', 1);
    plot(plotx, MC{ngroup}.p025, col2{ngroup}, 'LineWidth', 1);
    plot(plotx, MC{ngroup}.p975, col2{ngroup}, 'LineWidth', 1);
end
hold off;
    
title(titl);
xlabel(xlab);
ylabel(ylab);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
