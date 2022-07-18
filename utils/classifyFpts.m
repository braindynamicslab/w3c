function [stNodeIdx, unstNodeIdx, stSpiralIdx, sdlSpiralIdx, unstSpiralIdx] = classifyFpts(Jeig, tolr)
%CLASSIFYFPTS classify fixed points based on eigenvalues of their Jacobian
%matrices.
%   [stNodeIdx, unstNodeIdx, stSpiralIdx, unstSpiralIdx] = classifyFpts(Jeig)
% input:
%   Jeig: eigenvalues of Jacobian matries. N-by-d matrix, N is the number
%   of fixed points, d is the number of eigenvalues. 
%   tolr: tolerance. 0+-tolr is considered zero. 
% output:
%   stNodeIdx: indicator function of stable nodes (all eigenvalues are real 
%              and negative).
%   unstNodeIdx: indicator function of unstable nodes (all eigenvalues are
%              real and some are positive).
%   stSpiralIdx: indicator function of stable spirals (some eigenvalues
%               are complex, but all real parts are negative)
%   sdlSpiralIdx: indicator function of saddle spirals (the fixed point is
%               unstable in the subspace of real eigvectors, but stable in
%               its complement).
%   unstSpiralIdx: indicator function of unstable spirals (there is some
%               complex eigenvalue with a positive real part).
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 7-31-2019 ~
%}

% -- apply tolerance
if nargin>1 && ~isempty(tolr)
    Jeig = round(Jeig,floor(-log10(tolr)));
end
% -- classify fixed points
stNodeIdx=all(real(Jeig)<0,2) & all(imag(Jeig)==0,2);% stable nodes
unstNodeIdx=any(real(Jeig)>0,2) & all(imag(Jeig)==0,2);% unstable nodes
stSpiralIdx=all(real(Jeig)<0,2) & any(imag(Jeig)~=0,2);% stable spirals
sdlSpiralIdx=any(real(Jeig)>0 & imag(Jeig)==0,2) & any(real(Jeig)<0 & imag(Jeig)~=0,2);% saddle spirals
unstSpiralIdx=any(real(Jeig)>0 & imag(Jeig)~=0,2);% unstable spirals
end

