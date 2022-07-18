function [x,y]=mixmap(x1,y1,x2,y2)
% MIXMAP    fixing two maps y1=f1(x1), y2=f2(x2) to be y=f(x), where x
% contains elements from both x1 and x2, sorted in ascending order, such
% that y1=f(x1) and y2=f(x2). 
%{
~ Author: Mengsen Zhang <mengsenzhang@gmail.com> 07-17-2019 ~
%}
    x = [x1;x2];
    [x,idx] = sort(x);
    y = [y1;y2];
    y = y(idx,:);
end