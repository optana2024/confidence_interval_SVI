%check if the scalar x is included in the interval [c-r,c+r], up to a
%tolerance. If it is included, then t=1. Otherwise t=0
function t=isIncluded(x,c,r,tol)
if x >= c-r - tol && x <= c+r + tol
    t=1;
else
    t=0;
end
end
