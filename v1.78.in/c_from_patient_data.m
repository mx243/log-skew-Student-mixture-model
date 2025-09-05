function [c, namesout, averagevals, poststds, postmeans, prestds, premeans] ...
            = c_from_patient_data(patient_data, quadratic_mult, names, demeanfirst);
% Usage: [c, namesout, averagevals, poststds, postmeans, prestds, premeans] ...
%           = c_from_patient_data(patient_data, quadratic_mult, names, demeanfirst);
% This function takes in the raw data about the patients and outputs the
%  value of c (without the column of ones) to be used in MAIN.m.
% quadratic_mult should be a row vector of length Nrawvars,
%  containing a positive scalar if quadratic regression should be used for the
%  corresponding raw variable, and zero otherwise.
% Because (as Jensen says) the average of squares is greater than squares of averages 
%  (or equal if all contributors the same), a patient whose variables are all the average
%  value will not have all zeros after demeaning. averagevals is the vector of values
%  a patient with all variables average will have (so the squared ones will sometimes be negative).
% premeans is the set of means at the start, but will be zeros if there is no initial demeaning.
% prestds is the set of standard deviations before scaling, and will be ones if there is no initial descaling.
% premeans and prestds will have length corresponding to the unsquared variables only.
% postmeans is the set of means after adding the squared variables, and will be zeros if there
%  is no final demeaning.
% poststds is the set of stds after adding the squared variables, and will be ones if there
%  is no final descaling.
% postmeans and poststds will have length including squared variables also.
% If demeanfirst is 2, then:
%  demeaning and descaling will be done *both* before *and* after squaring;
%  averagevals will be zero on the unsquared variables and negative on the squared ones;
% If demeanfirst is 1, then: 
%  demeaning and destding will be done *before* squaring;
%  averagevals will be entirely zeros;
% if demeanfirst is 0, then::
%  squaring will be done before demeaning and descaling;
%  averagevals will be zero on the unsquared variables and negative on the squared ones.

% Change Log (as /home/rfs/ramakrishnan/software/lastact/tow24/c_from_patient_data.m):
%
%     1.1          25:sep:24    tow24    First version.
%     1.119        26:sep:24    rfs34    Made to work with old Matlab also, and outputs updated list of names.
%     1.130        07:oct:24    rfs34    averagevals and other outputs added.
%     1.133        08:oct:24    rfs34    demeanfirst option added, std now normalised by Npatients not Npatients-1.
%     1.137        10:oct:24    rfs34    demeanfirst=2 option added.
%     1.145        04:nov:24    rfs34    Added quadratic_mult option.
%     1.146        04:nov:24    rfs34    Added dithering for discrete variables.

% Change Log (as /home/rfs/ramakrishnan/software/lastact/mx243/c_from_patient_data.m):
%
%     1.67         04:nov:24    rfs34    Taken from tow24's code as an update.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.c_from_patient_data.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

use_quadratic = quadratic_mult > 0;

Npatients = size(patient_data, 1);
Nrawvars = size(patient_data, 2);
Nsquares = sum(use_quadratic);
Ntotvars = Nrawvars + Nsquares;

quadratic_mult = [ones(1, Nrawvars), quadratic_mult(quadratic_mult > 0)];

averagevals = mean(patient_data, 1);

if demeanfirst > 0,
   premeans = mean(patient_data, 1);
   prestds = std(patient_data, 1, 1);
else
   premeans = zeros(1, Nrawvars);
   prestds = ones(1, Nrawvars);
end

patient_data = patient_data - premeans(ones(Npatients, 1), :);
patient_data = patient_data ./ prestds(ones(Npatients, 1), :);
averagevals = averagevals - premeans;
averagevals = averagevals ./ prestds;

% Add the required squares.

patient_data_squared = patient_data .^ 2;

used_data_squared = patient_data_squared(:, use_quadratic);
used_data_squared = reshape(used_data_squared, [Npatients, Nsquares]);
patient_data = [patient_data, used_data_squared];

averagevalssquared = averagevals .^ 2;
used_averagevalssquared = averagevalssquared(:, use_quadratic);
averagevals = [averagevals, used_averagevalssquared];

namesout = names;
for nvar = 1 : Nrawvars;
   if use_quadratic(nvar),
      namesout = [namesout, {[names{nvar}, '^2']}];
   end
end

if demeanfirst ~= 1,
   postmeans = mean(patient_data, 1);
   poststds = std(patient_data, 1, 1);
else
   postmeans = zeros(1, Ntotvars);
   poststds = ones(1, Ntotvars);
end

patient_data = (patient_data - postmeans(ones(Npatients, 1), :)) ./ poststds(ones(Npatients, 1), :);
averagevals = (averagevals - postmeans) ./ poststds;

patient_data = patient_data .* quadratic_mult(ones(Npatients, 1), :);
averagevals = averagevals .* quadratic_mult;

% Add dithering for discrete variables.

for ntotvar = 1 : Ntotvars,
   v = patient_data(:, ntotvar);
   u = unique(v);
   if length(u) < 5,
      v = v + 0.01 * randn(size(v));
      patient_data(:, ntotvar) = v;
   end
end

c = patient_data;

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
