function OUT = local_celltomat(IN)
% convert cell-array to matrix

n =length(IN);

ok=1;
while ok
  if ~isempty(IN{ok})
    dims = size(IN{ok});
    ok = 0;
  else
    ok=ok+1;
  end
end

if length(dims)==2
  OUT = zeros(n,dims(1),dims(2));
  for k=1:n
    if ~isempty(IN{k})
      OUT(k,:,:) = IN{k};
    end
  end
  elseif length(dims)==3
  OUT = zeros(n,dims(1),dims(2),dims(3));
  for k=1:n
    if ~isempty(IN{k})
    OUT(k,:,:,:) = IN{k};
    end
  end
elseif length(dims)==4
  OUT = zeros(n,dims(1),dims(2),dims(3),dims(4));
  for k=1:n
    if ~isempty(IN{k})
      OUT(k,:,:,:,:) = IN{k};
    end
  end
elseif length(dims)==5
  OUT = zeros(n,dims(1),dims(2),dims(3),dims(4),dims(5));
  for k=1:n
    OUT(k,:,:,:,:,:) = IN{k};
  end
elseif length(dims)==6
  OUT = zeros(n,dims(1),dims(2),dims(3),dims(4),dims(5),dims(6));
  for k=1:n
    OUT(k,:,:,:,:,:,:) = IN{k};
  end
elseif length(dims)==7
  OUT = zeros(n,dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7));
  for k=1:n
    OUT(k,:,:,:,:,:,:,:) = IN{k};
  end
  
end

