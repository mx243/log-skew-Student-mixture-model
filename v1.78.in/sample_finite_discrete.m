function pos = sample_finite_discrete(log_p_vec, K)
% Usage: pos = sample_finite_discrete(log_p_vec, K)
% This function returns K samples from a finite discrete distribution. The
% row vec log_p_vec = log(p_vec) where p_vec is proportional to the vector
% of probabilities (some may be 0) of the different values the distribution 
% can take. 
% (We only consider the indices of these values, starting at 1; the vector of the values 
% themselves is available outside of this function.)

% Change Log:
%
%     1.1          26:aug:24    mx243      First version.
%     1.16         03:sep:24    rfs34      Sped up ~40-fold by vectorisation and minimising use of log domain.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_finite_discrete.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 2
    K = 1;
end

pos = find(log_p_vec == Inf, 1);
if ~isempty(pos),
   pos = pos(1, ones(1, K));
   return;
end

pos = NaN(1, K);

p_vec = exp(log_p_vec - max(log_p_vec(:)));

p_vec = p_vec(:);

% Having subtracted max(log_p_vec) there is no longer any need to be in the log domain.

sum_p_vec = sum(p_vec(:));

if sum_p_vec == 0,
   warning('sample_finite_discrete called with all zeros probabilities, returning all 1s');
   pos(:) = 1;
   return;
end

p_vec = p_vec / sum_p_vec;

cumsum_p_vec = cumsum(p_vec);

u = rand(1, K);

pos = 1 + sum(u(ones(length(p_vec), 1), :) > cumsum_p_vec(:, ones(1, K)), 1);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
