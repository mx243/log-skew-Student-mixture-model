function plot_bar(ax, values, val_true)
% Usage: plot_bar(ax, values, val_true)
% This function plots the histogram of a set of values and the true value
% on the figure specified by ax.

% Change Log:
%
%     1.1          30:jul:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_bar.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

values = squeeze(values);
if size(values, 1) > 1,
    values = values';
end

left = min(val_true, min(values));
right = max(val_true, max(values));
tmpx = left : (right - left) / 50 : right;

if(isempty(tmpx))
    tmpx = [left, left + 1];
end

tmpy = hist(values, tmpx);
nvtrue = hist(val_true * ones(1, max(tmpy)), tmpx);

h = bar(ax, tmpx, [tmpy(:), nvtrue(:)]);

set(h(1), 'FaceColor', [0, 0, 1]);
set(h(2), 'FaceColor', [0, 1, 0]);

set(h(1), 'EdgeColor', 'none');
set(h(2), 'EdgeColor', 'none');

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
