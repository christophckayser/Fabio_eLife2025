
function out = local_combinebehavpart(Behav)

% for each field we cast participants into the first dimension
% and combine data
fieldnames = fields(Behav(1));
nsub = length(Behav);
out = [];
for f=1:length(fieldnames)
  tmp = getfield(Behav,{1},fieldnames{f});
  tmp = zeros([1 size(tmp)]);
  for S=1:nsub
    tmp(S,:,:) = getfield(Behav,{S},fieldnames{f});
  end
  out  = setfield(out,fieldnames{f},tmp);
end