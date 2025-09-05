function plotphi(samplefile);
% Usage: plotphi(samplefile);
% This function plots the trajectories of elements of phi without reference 
%  to what is adherent to each cluster. One should therefore expect the 
%  various trajectories to converge not necessarily to their own true value,
%  but to one of the true values, and the colours to correspond between plots.

% Change Log:
%
%     1.1          18:sep:24    rfs      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plotphi.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

load(samplefile);

N = size(Samples.nonempty_c_num, 2);
K = size(Samples.alpha, 1);
D = Priors.D;
smpl_cnt_vec = 1 : N;

[psamples{1 : N}] = deal(Samples.phi.p);
psam = arrangesamples(psamples);
plotvec(psam, Args_true.phi.p);
title('p');

[msamples{1 : N}] = deal(Samples.phi.m);
msam = arrangesamples(msamples);
plotvec(msam, Args_true.phi.m);
title('m');

[musamples{1 : N}] = deal(Samples.phi.mu);
musam = arrangevectorsamples(musamples);
for d = 1 : D,
   plotvec(reshape(musam(d, :, :), [size(musam, 2), size(musam, 3)]), Args_true.phi.mu(d, :));
   title(sprintf('$\\mu_{%d}$', d), 'Interpreter', 'Latex');
end

[Ssamples{1 : N}] = deal(Samples.phi.S);
Ssam = arrangematrixsamples(Ssamples);
for d1 = 1 : min(D, 3),
   for d2 = d1 : min(D, 3),
      plotvec(reshape(Ssam(d1, d2, :, :), [size(Ssam, 3), size(Ssam, 4)]), squeeze(Args_true.phi.S(d1, d2, :)).');
      title(sprintf('$S_{%d, %d}$', d1, d2), 'Interpreter', 'Latex');
   end
end

[nusamples{1 : N}] = deal(Samples.phi.nu);
nusam = arrangevectorsamples(nusamples);
for d = 1 : D,
   plotvec(reshape(nusam(d, :, :), [size(nusam, 2), size(nusam, 3)]), Args_true.phi.nu(d, :));
   title(sprintf('$\\nu_{%d}$', d), 'Interpreter', 'Latex');
end

return;


function [out] = arrangematrixsamples(samples);
% Usage: [out] = arrangematrixsamples(samples);
% samples is a cell array with each element being a 3-d of various size(..., 3) values.
% out is then a rectanguloid array with NaN where the arrays were too short in dimension 3.

D = size(samples{1}, 1);
sz3 = @(c) size(c, 3);
C = max(cellfun(sz3, samples));
N = length(samples);
out = NaN(D, D, C, N);
for n = 1 : N,
   l = sz3(samples{n});
   out(:, :, 1 : l, n) = samples{n};
end

return;


function [out] = arrangevectorsamples(samples);
% Usage: [out] = arrangevectorsamples(samples);
% samples is a cell array with each element being a 2-d of various size(..., 2) values.
% out is then a rectanguloid array with NaN where the arrays were too narrow.

D = size(samples{1}, 1);
sz2 = @(c) size(c, 2);
C = max(cellfun(sz2, samples));
N = length(samples);
out = NaN(D, C, N);
for n = 1 : N,
   l = sz2(samples{n});
   out(:, 1 : l, n) = samples{n};
end

return;


function [out] = arrangesamples(samples);
% Usage: [out] = arrangesamples(samples);
% samples is a cell array with each element being a column vector of various lengths.
% out is then a rectangular array with NaN where the column vectors were too short.

C = max(cellfun(@length, samples));
N = length(samples);
out = NaN(C, N);
for n = 1 : N,
   l = length(samples{n});
   out(1 : l, n) = samples{n};
end

return;


function plotvec(samples, truth);
% Usage: plotvec(samples, truth);
% Used for plotting each row of samples in a different colour,
%  corresponding to each row of truth.

cols = 'rgbcmyk';
Ncols = length(cols);

truth = [truth, NaN(1, size(samples, 1) - length(truth))];

figure;
clf;
hold on;
for k = 1 : min(Ncols, size(samples, 1)),
   plot([1 : size(samples, 2)], samples(k, :), [cols(k), '.']);
end
for k = 1 : min(Ncols, size(samples, 1)),
   plot([1, size(samples, 2)], truth([k, k]), [cols(k), '-x'], 'LineWidth', 1.5, 'MarkerSize', 10);
end
xlabel('Sample number');
hold off;

return;




% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
