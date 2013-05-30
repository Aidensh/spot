function x = lsqrdivide(A,b,mode)

if mode == 2
    A = A';
end

[m,n] = size(A);
opts = spotparams;
maxits = opts.cgitsfact * min(m,min(n,20));

% Preallocate x
x(size(A,2),size(b,2)) = cast(0,class(b));
for j=1:size(b,2)
    x(:,j) = spot.solvers.lsqr(m,n,A,b(:,j),opts.cgdamp,opts.cgtol,...
        opts.cgtol,opts.conlim,maxits,opts.cgshow);
end