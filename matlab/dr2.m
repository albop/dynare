% Copyright (C) 2001 Michel Juillard
%
% function used for dr_algo == 1
function ghs2=dr2(ys)
  global fname_ dr_
  
  dr_.ys = ys;
  fh = str2func([fname_ '_fff']);
  dr_.fbias = 2*feval(fh,dr.ys);
  dr1(0);
  ghs2 = dr_.ghs2;
