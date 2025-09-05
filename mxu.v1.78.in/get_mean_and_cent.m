function [MClifetime, MCsurvival, MChazardrt] = get_mean_and_cent(lifetime, survival, hazardrt, usealpha)
% Usage: [MClifetime, MCsurvival, MChazardrt] = get_mean_and_cent(lifetime, survival, hazardrt, usealpha)
% This function returns the means and centiles of the curves lifetime,
% survival, hazardrt.
%
% In the case usealpha ~= 0, the input lifetime will actually be the
% survival data, and survival and hazardrt are meaningless (in fact empty).
% Hence MCsurvival and MChazardrt are set to [].

% Change Log:
%
%     1.1          10:oct:24    mx243      First version.
%     1.56         17:oct:24    mx243      Added usealpha.
%     1.57         22:oct:24    mx243      Switched to using ceil(0.975 * Ncurves) for the 97.5 centile
%                                          instead of floor.
%     1.71         25:jan:25    rfs34      Comments only changed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.get_mean_and_cent.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

Ncurves = size(lifetime, 1);
Nx = size(lifetime, 2);

MClifetime = struct('mean', NaN(1, Nx), 'p025', NaN(1, Nx), 'p975', NaN(1, Nx));
if ~usealpha
    MCsurvival = struct('mean', NaN(1, Nx), 'p025', NaN(1, Nx), 'p975', NaN(1, Nx));
    MChazardrt = struct('mean', NaN(1, Nx), 'p025', NaN(1, Nx), 'p975', NaN(1, Nx));
else
    MCsurvival = [];
    MChazardrt = [];
end

LB = ceil(0.025 * Ncurves);
UB = ceil(0.975 * Ncurves);

lifetime = sort(lifetime, 1);
MClifetime.mean = mean(lifetime, 1);
MClifetime.p025 = lifetime(LB, :);
MClifetime.p975 = lifetime(UB, :);

if ~usealpha
    survival = sort(survival, 1);
    MCsurvival.mean = mean(survival, 1);
    MCsurvival.p025 = survival(LB, :);
    MCsurvival.p975 = survival(UB, :);

    hazardrt = sort(hazardrt, 1);
    MChazardrt.mean = mean(hazardrt, 1);
    MChazardrt.p025 = hazardrt(LB, :);
    MChazardrt.p975 = hazardrt(UB, :);
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
