function [dnx,dny]=grad(n,x,y)

dnx=(n(2:end-1,3:end)-n(2:end-1,1:end-2))./(x(2:end-1,3:end)-x(2:end-1,1:end-2));
dny=(n(3:end,2:end-1)-n(1:end-2,2:end-1))./(y(3:end,2:end-1)-y(1:end-2,2:end-1));


