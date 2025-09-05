function y = piecewise_exponential(points, vals, grads, bounds)
% Usage: y = piecewise_exponential(points, vals, grads, bounds)
% Samples from piecewise exponential (exponential components defined by
%  the equations of their log). 
% bounds: array of endpoints of the components.
% points, vals, grads: (points(j), vals(j)) lie on the logarithm of the jth
%  component, which has gradient grads(j).
% y has size [1, 2]; y(1) is the desired sample, y(2) is component number it came from.

% Change Log:
%
%     1.1          28:sep:19    jc2062   As first received from JC2062.
%     1.2          10:oct:19    jc2062   As received from jc2062.
%     1.10         12:oct:19    rfs34    Dimension passed to log_total.
%     1.17.3.3     15:oct:19    rfs34    Branch used with prior like gingernut, big Nsamples, 
%                                        and for debugging this ars.m and piecewise_exponential.m .
%     1.18         15:oct:19    rfs34    This branch of this function and ars.m taken back into main line.
%     1.21         16:oct:19    rfs34    Now catches underflow due to very negative vals.
%     1.28.1.2     22:oct:19    rfs34    Attempts to catch problems at the end.
%     1.28.1.3     22:oct:19    rfs34    Now attempts to catch infinite grads.
%     1.28.1.4     22:oct:19    rfs34    Trivial bug in above fixed.
%     1.28.1.5     22:oct:19    rfs34    Catches infinite grads earlier; catches NaN weights
%     1.39         22:oct:19    rfs34    Taken this code back into the main line.
%     1.43         23:oct:19    rfs34    Comments only changed.

% Change Log (as ~rfs/ramakrishnan/software/lastact/mx243/piecewise_exponential.m):
%
%     1.1          08:jul:24    mx243    Taken from above, and modified by rfs34 to avoid
%                                        using i as an iterator and improve layout.
%     1.69         06:nov:24    rfs34    Modified to prevent multiple infinite weights leading
%                                        to always picking the top segment.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.piecewise_exponential.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

maxvals = max(vals);
if isinf(maxvals),
   if maxvals == -Inf,
      warning('All exp(vals) are zero in piecewise_exponential');
   else
      error('Infinite positive vals in piecewise_exponential');
   end
end
vals = vals - maxvals;

shifted_bounds = bounds(2 : end);
bounds = bounds(1 : end - 1);

maxbordervals = max([grads .* (bounds - points) + vals, ...
                     grads(end) .* (shifted_bounds(end) - points(end)) + vals(end), ...
                     0]);
vals = vals - maxbordervals;

index = find(grads);

weights = zeros([1, length(bounds)]);
weights(index) = (1 ./ grads(index)) ...
                 .* (exp(grads(index) .* (shifted_bounds(index) - points(index)) + vals(index)) ...
                     - exp(grads(index) .* (bounds(index) - points(index)) + vals(index)));
zindex = find(~grads);
weights(zindex) = exp(vals(zindex)) .* (shifted_bounds(zindex) - bounds(zindex));

% RFS: Added this because we were getting NaN weights getting picked.
weights(isnan(weights)) = 0;

z = general_discrete(weights);
u = rand;

bz = bounds(z);
sbz = shifted_bounds(z);

if grads(z) == Inf,
   y = [sbz, z];
elseif grads(z) == -Inf,
   y = [bz, z];
elseif isinf(bz) && grads(z) <= 0,
   % RFS: This shouldn't be able to happen, and should have been checked for in calling function.
   warning('Piecewise exponential called with unbounded density at left end');
   y = [sbz * (1 - 0.1 * sign(sbz)), z];
elseif isinf(sbz) && grads(z) >= 0,
   warning('Piecewise exponential called with unbounded density at right end');
   y = [bz * (1 + 0.1 * sign(bz)), z];
elseif grads(z) == 0,
   y = [bz .* u + sbz .* (1 - u), z];
else
    a = [log(u) + (grads(z) * sbz), log(1 - u) + grads(z) * bz];
    y = [(1 / grads(z)) * log_total(a, 2), z];
end

% RFS: I added this because otherwise numerical inaccuracy in the last two lines above
%      was leading to these conditions being violated.
if y(1) <= bz,
   y(1) = min([0.001 * sbz + 0.999 * bz, ...
               bz * (1 + sign(bz) * 0.001)]);
end
if y(1) >= sbz,
   y(1) = max([0.001 * bz + 0.999 * sbz, ...
              sbz * (1 - sign(sbz) * 0.001)]);
end

if any(isinf(y)) || any(isnan(y)) || y(1) < bounds(y(2)) || y(1) > shifted_bounds(y(2)),
   keyboard;
end

end

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
