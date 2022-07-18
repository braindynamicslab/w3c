function v = belowdiag(X)
%BELOWDIAG Elements below the diagonal of a square matrix X linearized to a
%vector v. 
v = X(tril(true(size(X)),-1));
end

