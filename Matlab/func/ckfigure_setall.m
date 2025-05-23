function ckfigure_setall(handle,prop,arg)

% script to set a property for all subplots of a given figure
%  ckfigure_setall(handle,prop,arg)
%

Ch = get(handle,'Children');

for l=1:length(Ch)
  try
    set(Ch(l),prop,arg);
  catch 
  end
end
