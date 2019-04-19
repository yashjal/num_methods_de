function u = be(A,un,f,dt)
% backward euler step in time
% f = f(xn+1,tn+1)
I = speye(size(A));
X = (I-dt*A);
b = un+dt*f;

[R,flag] = chol(X);
if ~flag
    u = R\(R'\b);
else
    disp('Factorization failed.');
    u = -X\b;
end

end