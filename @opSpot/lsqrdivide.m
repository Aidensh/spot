function x = lsqrdivide(A,b)

[m,n] = size(A);
opts = spotparams;
maxits = opts.cgitsfact * min(m,min(n,20));

% Preallocate x
x = zeros(size(b));
for j=1:size(b,2)
    x(:,j) = spot.solvers.lsqr(m,n,A,b(:,j),opts.cgdamp,opts.cgtol,...
        opts.cgtol,opts.conlim,maxits,opts.cgshow);
end