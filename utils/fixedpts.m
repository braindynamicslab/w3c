function [fpts,Jeig,flag,idxICs] = fixedpts(fhandle,Xs,tolr,opt)
%FIXEDPTSGRID find fixed points given a set of initial conditions
%   [fpts,Jeig,flag] = fixedpts(fhandle,Xs,tolr)
% input:
%   fhandle: handle of a function the zeros of which we seek
%   Xs: N-by-d matrix, specifying initial conditions for zero-search. N is
%   the number of initial conditions, and d is the dimension of the domain
%   of the function. 
%   tolr: a number that defines how close two fixed points are to be
%   considered as identical. 
% output:
%   fpts: coordinates of fixed points, 1 pt = 1 row Jeig: eigenvalue of the
%   Jacobian matrix at the fixed points. 
%   flag: the flag returned by fsolve (note that all returns flagged as
%   "no-solution" has been disgarded.
%   idx: index of the initial condition that leads to the solution. 
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com 7-15-2019 ~
%}

if nargin<4
    opt = optimset('Display','none');
end

% - tolerance (for treating two solutions as the same)
if nargin<3
    tolr=1e-13;
end

[Nx0,d]=size(Xs);

[fpts,Jeig]=deal(nan(Nx0,d));
flag=nan(Nx0,1);


% -- find zeros of 2d function
for n=1:Nx0
    [fpts(n,:),~,flag(n),~,J]=fsolve(fhandle,Xs(n,:)',opt);
    Jeig(n,:)=eig(J);
    clear J
end

fpts=fpts(flag>0,:);
Jeig=Jeig(flag>0,:);
[~,idx]=unique(round(fpts,floor(-log10(tolr))),'rows');

fpts=fpts(idx,:);
Jeig=Jeig(idx,:);
idxICs=find(flag>0);
idxICs=idxICs(idx);
flag=flag(idxICs);

end

