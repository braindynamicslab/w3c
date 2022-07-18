function bwX = bwValues(X)
%BWVALUES values between rows of a matrix.
%   bwX = bwValues(X)
% input:
%   X: N-by-M matrix.
% output:
%   bwX: (N-1)-by-M matrix, where the n-th row is the average of the n-th
%   and the (n+1)-th row of the input matrix.
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 7-15-2019 ~
%}
bwX = (X(1:end-1,:) + X(2:end,:))/2; 
end

