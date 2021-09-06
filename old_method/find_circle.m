function [nc,r0,the1]=find_circle(N,x1,y1,x0,y0,x_B0,y_B0,R0,the0,nc,vv,dt)

r0=zeros(1,N);
the1=the0;
for ii=1:N
    index=find((x1(ii)-x_B0).^2+(y1(ii)-y_B0).^2<=R0.^2);
    if ~isempty(index)
        nc(ii)=index;
        r0(ii)=sqrt((x1(ii)-x_B0(index)).^2+(y1(ii)-y_B0(index)).^2);
        the1(ii)=atan2(y1(ii)-y_B0(index),x1(ii)-x_B0(index));
    else
        if nc(ii)~=0
            r0(ii)=sqrt((x0(ii)-x_B0(nc(ii))).^2+(y0(ii)-y_B0(nc(ii))).^2);%粒子在圈内（圆周运动后、碰撞到圈外前）相对圆心的角度
            ww=vv(ii)./r0(ii);
            dthe=ww*dt;
            the1(ii)=the0(ii)+dthe;
            
        end 
        nc(ii)=0;
    end 
    
end