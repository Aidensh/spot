function x = lsqrdivide(A,b,mode)

if mode == 2
    A = A';
end

[m,n]  = size(A);
opts   = spotparams;
maxits = opts.cgitsfact * min(m,min(n,20));

% Run lsqr solver
for j=size(b,2):-1:1
  if opts.uselsmr % run lsmr instead of lsqr
    x(:,j) = spot.solvers.lsmr(A,b(:,j),opts.cgdamp,opts.cgtol,...
        opts.cgtol,opts.conlim,maxits,opts.reorthosize,opts.cgshow);
  else
    x(:,j) = spot.solvers.lsqr(m,n,A,b(:,j),opts.cgdamp,opts.cgtol,...
        opts.cgtol,opts.conlim,maxits,opts.cgshow);
  end
end