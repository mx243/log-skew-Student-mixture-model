function ctrl_file = write_Stu_mix_ctrl(save_dir, save_id, reverse, speedup)
% Usage: ctrl_file = write_Stu_mix_ctrl(save_dir, save_id, reverse, speedup)
% This function creates a .mat file specifying a set of ctrl arguments for
% SA in the Student mixture model.
% If reverse is passed and is 1 then the schedule is reversed.
% If it is passed and is 2 then the reversed schedule is tacked on to the end 
%  of the forward schedule.
% If speedup is passed then annealing is sped up by increasing the gaps between
%  coolness values by a factor of speedup.

% Change Log:
%
%     1.1          04:sep:24    mx243      First version.
%     1.20         09:sep:24    rfs34      Number of annealing samples reduced to 1000.
%     1.22         10:sep:24    mx243      Added return value ctrl_file, added -v6 to the save command.
%     1.29         16:sep:24    rfs34      Reduced annealing time while debugging.
%     1.31         16:sep:24    rfs34      Back to 1000 annealing steps.
%     1.33         17:sep:24    rfs34      Added a run-in at zero coolness to the start.
%     1.36         18:sep:24    rfs34      Changed annealing schedule to exponential.
%     1.37         20:sep:24    rfs34      Added discards and coarsened schedule in compensation.
%     1.38         20:sep:24    rfs34      Removed the coarsening.
%     1.36         20:sep:24    rfs34      Reverted to settings of v1.36 as the above didn't help.
%     1.41         23:sep:24    rfs34      Switched to a finer schedule with most fineness on the high coolness end;
%                                          added reversal option; added 2 discards for every sample.
%     1.42         24:sep:24    rfs34      Fixed typo in above.
%     1.43         24:sep:24    rfs34      Removed discards; added reverse=2 option.
%     1.48         03:oct:24    rfs34      Added speedup option.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.write_Stu_mix_ctrl.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if nargin < 3 || isempty(reverse),
   reverse = 0;
end
if nargin < 4 || isempty(speedup),
   speedup = 1;
end

% t_SA = [0, exp([-10 : 0.01 : 0])];

logoddst = [[-15 : (0.04 * speedup) : 0], [0.01 : (0.01 * speedup) : 10]];
t_SA = [0, exp(logoddst) ./ (1 + exp(logoddst)), 1, 1];

% Point after which samples at coolness 1 can be inserted.
Nschedinsert = length(t_SA);

if reverse == 1,
   t_SA = flipdim(t_SA, 2);
elseif reverse == 2,
   t_SA = [t_SA, flipdim(t_SA, 2)];
end

% If reverse == 2 then we average over annealing in the two directions
% to try to remove the effect of convergence lag during annealing.
dirsign = sign([0, diff(t_SA)]) / (1 + (reverse == 2));

t_SA_num = length(t_SA);
nsmpl = [30, 1 * ones(1, t_SA_num - 1)];
dis = [29, 0 * ones(1, t_SA_num - 1)];

ctrl_file = sprintf('%s/sa_ctrl_Stu_mix.%s.mat', save_dir, save_id);
save(sprintf('%s/sa_ctrl_Stu_mix.%s.mat', save_dir, save_id), '-v6');

return;

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
