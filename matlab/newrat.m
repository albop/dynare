%
%function [xparam1, x, h, g, fval0] = newrat(func0,x,hh,gg,varargin)
function [xparam1, hh, gg, fval] = newrat(func0,x,hh,gg,varargin)


icount=0;
nx=length(x);
xparam1=x;
lamtol=1.e-7;
ftol=1.e-5;
options=optimset('fminunc');
options.MaxFunEvals=200;
options.TolFun= 1.0000e-005;
options.MaxIter=1;
options.LargeScale='off';
func = str2func(func0);
fval0=feval(func,x,varargin{:});
if isempty(hh)
    [dum, gg, t3, s3, x3, h1]=mr_hessian(func0,x,varargin{:});
    hh = reshape(dum,nx,nx);
end
disp(['Gradient norm ',num2str(norm(gg))])
disp(['Minimum Hessian eigenvalue ',num2str(min(eig(hh)))])
disp(['Maximum Hessian eigenvalue ',num2str(max(eig(hh)))])
g=gg;
h{1}=hh;
check=0;
if max(eig(hh))<0, disp('Negative definite Hessian! Local maximum!'), pause, end,
while norm(gg)>1.e-3 & check==0,
    icount=icount+1;
    disp([' '])
    disp(['Iteration ',num2str(icount)])
    x0=xparam1-inv(hh)*gg;
    fval=feval(func,x0,varargin{:});
    c=mr_nlincon(x0,varargin{:},1);
    lam=1;
    while c
        lam=lam*0.9;
        x0=xparam1-inv(hh)*gg.*lam;
        fval=feval(func,x0,varargin{:});
        c=mr_nlincon(x0,varargin{:},1);
    end        
    if (fval0(icount)-fval)<ftol,
        disp('Try line search')
        [lam,fval,EXITFLAG,OUTPUT,GRAD,HESSIAN]=fminunc(@lsearch, 0, options, func, xparam1, inv(hh)*gg , varargin{:});
        x0=xparam1-inv(hh)*gg.*lam;
%         fvala=fval;
%         x0a=x0;
    end
%     if (fval0(icount)-fval)<ftol & min(eig(hh))<0,
%         disp('Try direction of largest negative eigenvalue')
%         [v, dum]=eig(hh);
%         [id, ij]=min(diag(dum));
%         [lam,fval,EXITFLAG,OUTPUT,GRAD,HESSIAN]=fminunc(@lsearch, 0, options, func, xparam1, v(:,ij) , varargin{:});
%         x0=xparam1-v(:,ij)*lam;
%     end
    if (fval0(icount)-fval)<ftol*ftol,
%         if fvala<fval,
%             fval=fvala;
%             x0=x0a;
%         end
        disp('No further improvement is possible!')
        check=1;
    end
    
    
    xparam1=x0;
    x(:,icount+1)=xparam1;
    fval0(icount+1)=fval;
    disp(['LAMBDA        ',num2str(lam)])
    disp(['DX norm       ',num2str(norm(inv(hh)*gg.*lam))])
    disp(['FVAL          ',num2str(fval)])
    disp(['Improvement   ',num2str(fval0(icount)-fval)])
    
    if norm(x(:,icount)-xparam1)>1.e-12,
        %[dum, gg, t3, s3, x3, h1]=hessian('mj_optmumlik',xparam1,gend,data,1);
        [dum, gg, t3, s3, x3, h1]=mr_hessian(func0,xparam1,varargin{:});
        hh = reshape(dum,nx,nx);
    end
    disp(['Gradient norm  ',num2str(norm(gg))])
    disp(['Minimum Hessian eigenvalue ',num2str(min(eig(hh)))])
    disp(['Maximum Hessian eigenvalue ',num2str(max(eig(hh)))])
    if max(eig(hh))<0, disp('Negative definite Hessian! Local maximum!'), pause, end,
    
    h{icount+1}=hh;
    g(:,icount+1)=gg;
    save m1 x h g fval0
end

return

%  
function f00 = lsearch(lam,func,x,dx,varargin)


x0=x-dx*lam;
f00=feval(func,x0,varargin{:});





