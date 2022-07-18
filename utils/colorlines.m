function colorlines(lh,c)
%COLORLINES color lines given the handle of lines "lh", and colors "c" (a
%N-by-3 matrix with N = length(lh);
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 07-03-2019 ~
%}
N=length(lh);
cmap=mat2cell(c,ones(N,1),3);
[lh.Color]=deal(cmap{:});
end

