function [A,B,ys] = dynare_resolve()
  global dr1_test_ ys_ dr_
  
  dr_ = resol1(ys_,0,1,1);

  if dr1_test_(1) > 0
    A = [];
    B = [];
    ys = [];
    return
  end
  
  [A,B] = kalman_transition_matrix(dr_);
  ys = dr_.ys;
