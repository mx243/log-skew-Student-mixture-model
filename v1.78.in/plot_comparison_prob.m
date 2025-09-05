function plot_comparison_prob(plotx, ploty, Ngroup, titl)
% Usage: plot_comparison_prob(plotx, ploty, Ngroup, titl)
% This function plots the mean of the posterior distribution of the
% probability that group 1 has larger ploty than group 2, with the uniform
% distr on [0, 1] as prior.

% Change Log:
%
%     1.1          23:oct:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_comparison_prob.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if ~(Ngroup == 2 && size(ploty, 2) == 2)
    error('plot_comparison_prob called with one set of data.');
end

Nbigger = zeros(size(plotx));
Nequal = zeros(size(plotx));

Nrows1 = size(ploty{1}, 1);
Nrows2 = size(ploty{2}, 1);

Npairs = Nrows1 * Nrows2;

for it_N1 = 1 : Nrows1
    for it_N2 = 1 : Nrows2
        Nbigger = Nbigger + (ploty{1}(it_N1, :) > ploty{2}(it_N2, :));
        Nequal = Nequal + (ploty{1}(it_N1, :) == ploty{2}(it_N2, :));
    end
end

probs = (Nbigger + Nequal / 2 + 1) / (Npairs + 2); % The mean of said posterior at the beginning.

figure;
plot(plotx, probs, 'b-');

title(['Probability of group 1 having higher ', titl, ' than group 2 against time']);
xlabel('Time(days)');
ylabel('Probability');
ylim([-0.01, 1.01]);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
