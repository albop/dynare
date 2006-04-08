% solves x-a*x*a'=b for b (and then x) symmetrical
function [x,ns_var]=lyapunov_symm(a,b)
  global options_ 
  
  info = 0;
  if size(a,1) == 1
    x=b/(1-a*a);
    return
  end
  [u,t] = schur(a);
  e1 = abs(ordeig(t)>2-options_.qz_criterium);
  k = sum(e1);
  [u,t] = ordschur(u,t,e1); 
  n = length(e1)-k;
  b=u(:,k+1:end)'*b*u(:,k+1:end);
  t = t(k+1:end,k+1:end);
  x=zeros(n,n);
  for i=n:-1:2
    if t(i,i-1) == 0
      if i == n
	c = zeros(n,1);
      else
	c = t(1:i,:)*(x(:,i+1:end)*t(i,i+1:end)')+...
	    t(i,i)*t(1:i,i+1:end)*x(i+1:end,i);
      end
      q = eye(i)-t(1:i,1:i)*t(i,i);
      x(1:i,i) = q\(b(1:i,i)+c);
      x(i,1:i-1) = x(1:i-1,i)';
    else
      if i == n
	c = zeros(n,1);
	c1 = zeros(n,1);
      else
	c = t(1:i,:)*(x(:,i+1:end)*t(i,i+1:end)')+...
	    t(i,i)*t(1:i,i+1:end)*x(i+1:end,i)+...
	    t(i,i-1)*t(1:i,i+1:end)*x(i+1:end,i-1);
	c1 = t(1:i,:)*(x(:,i+1:end)*t(i-1,i+1:end)')+...
	     t(i-1,i-1)*t(1:i,i+1:end)*x(i+1:end,i-1)+...
	     t(i-1,i)*t(1:i,i+1:end)*x(i+1:end,i);
      end
      q = [eye(i)-t(1:i,1:i)*t(i,i) -t(1:i,1:i)*t(i,i-1);...
	   -t(1:i,1:i)*t(i-1,i) eye(i)-t(1:i,1:i)*t(i-1,i-1)];
      z =  q\[b(1:i,i)+c;b(1:i,i-1)+c1];
      x(1:i,i) = z(1:i);
      x(1:i,i-1) = z(i+1:end);
      x(i,1:i-1)=x(1:i-1,i)';
      x(i-1,1:i-2)=x(1:i-2,i-1)';
      i = i - 1;
    end
  end
  if i == 2
    c = t(1,:)*(x(:,2:end)*t(1,2:end)')+t(1,1)*t(1,2:end)*x(2:end,1);
    x(1,1)=(b(1,1)+c)/(1-t(1,1)*t(1,1));
  end
  x=u(:,k+1:end)*x*u(:,k+1:end)';
  ns_var = [];
  ns_var = find(any(abs(u(:,1:k)) > 1e-8,2)); 
