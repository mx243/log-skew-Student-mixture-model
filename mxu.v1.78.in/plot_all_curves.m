function plot_all_curves(plotx, ploty, Ngroups, ord, color, titl, xlab, ylab)
% Usage: plot_all_curves(plotx, ploty, Ngroups, ord, color, titl, xlab, ylab)
% This function plots (possibly more then one set of) samples in ploty. 
% The ord it takes in randomises the order in which the samples are plotted. 
% Precisely:
% plotx should be a vector of size [1, Nx].
% ploty should be have size {1, Ngroups};
%  ploty{ngroup} should have size [Ncurves(ngroup), Nx].
% ord should be a vector of length sum(Ncurves), containing 
%  for each ngroup, Ncurves(ngroup) elements of value ngroup.
%  The order in which these occur controls the order in which 
%  the curves of the various groups are plotted.
% color should have size {1, Ngroups}, with color{ngroup} being
%  the colour and line style (e.g. 'r--' for a red dashed line)
%  in which curves of that group are to be plotted.
% titl, xlab, and ylab are respectively the desired title, xlabel, and ylabel for the plot.

% Change Log:
%
%     1.1          10:oct:24    mx243      First version.
%     1.52         14:oct:24    rfs34      Comments extended.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_all_curves.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

pos = zeros(1, Ngroups); % Indicates which row to plot for each set of samples.

figure;
for num = 1 : size(ord, 2)
    hold on;
    pos(ord(num)) = pos(ord(num)) + 1;
    plot(plotx, ploty{ord(num)}(pos(ord(num)), :), color{ord(num)}); 
end

title(titl);
xlabel(xlab);
ylabel(ylab);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
