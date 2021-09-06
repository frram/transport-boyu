%d: distance between two cell centers
%D0: pure diffusion coefficient
%v_p: v_parallel
%v_r: v_random
%n: n=d/(v_p*dt) indicates step length of parallel velocity
%m: %v_r/v_p=1/mm indicates step length of random velocity
%Pe: Peclet number. Pe = v_p*d/D0 = v_p*d/(v_r^2*dt) =
%(v_p/v_r)^2*d/(v_p*dt) = m^2*n
%rat=R0/width indicates the gap width

%nc: the number of cells in which the particles are located. if nc=0, then
%out of cell.
%r0 & the0: polar coordinates of particles in the cell
%x0 & y0: vertical coordinates of partilces for ploting figure

%% fix d,D0,m，(i.e., fix v_p/v_r, scan step length)
clear;clc
d=1;
D0=0.1;
m=5;
tau0=d^2/D0;

%scan radio
rat=[10:10:150];
kk_r=0;
for rrat=rat
    kk_r=kk_r+1;
    
    %scan cell numbers for every line and column
    num=3:2:15;
    kk_num=0;
    for nnum=num
        kk_num=kk_num+1;
        l=5*d/nnum;
        RR=l/2*rrat./(rrat+1);
        wwid=l/2*1./(rrat+1);
        
        % set magnetic cells: centers and radii
        x_B0=-(nnum-1)/2*l:l:(nnum-1)/2*l;
        y_B0=-(nnum-1)/2*l:l:(nnum-1)/2*l;
        [x_B0,y_B0]=meshgrid(x_B0,y_B0);
        x_B0=x_B0(:)';
        y_B0=y_B0(:)';
        n0=length(x_B0);
        nc_boundry=find(x_B0==min(x_B0) | x_B0 == max(x_B0) | y_B0==min(y_B0) | y_B0 == max(y_B0));
        RR0=ones(1,n0)*RR;
        
        %set phi and magnetic field
        x=-2.5*d:0.01*d:2.5*d;
        y=-2.5*d:0.01*d:2.5*d;
        [x,y]=meshgrid(x,y);
        phi=zeros(size(x));
        for ii=1:n0
            index=find((x-x_B0(ii)).^2+(y-y_B0(ii)).^2<=RR0(ii)^2);
            phi(index)=cos(pi*(x(index)-x_B0(ii))/RR0(ii)/2).*cos(pi*(y(index)-y_B0(ii))/RR0(ii)/2);
        end
        [By,Bx]=grad(phi,x,y);
        x_net=x(2:end-1,2:end-1);
        y_net=y(2:end-1,2:end-1);
        By=By;
        Bx=-Bx;
        
        %scan n (Pe)
        kk_n=0;
        n=10;
        for nn=n
            kk_n=kk_n+1;
            clear t_count;
            clear N1;
            
            PPe=m^2*nn;
            v_p=m^2*nn*D0/d;
            v_r=m*nn*D0/d;
            dt=tau0/nn/nn/m/m;
            dl=v_r*dt;
            
            tau_eff=tau0/sqrt(PPe);
            t=0:dt:6*tau_eff;
            
            clear D_eff_tem
            for ll=1:3 %averaging
                %initialize the particle coordinates
                N=2000;
                nc=(n0+1)/2*ones(1,N);
                r0=rand(1,N).*RR0(nc);
                the0=rand(1,N)*2*pi-pi;
                x0=r0.*cos(the0)+x_B0(nc);
                y0=r0.*sin(the0)+y_B0(nc);
                
                v_net=-5*v_p:0.01*v_p:5*v_p;%set parallel velocity distribution
                fv=1/pi/v_p*exp(-v_net.^2/pi/v_p^2);
                vv=dist(1,N,v_net,fv);
                % vv=abs(vv);
                
                % %check the velocity distribution
                % vx0=vv.*sin(the0);
                % vy0=-vv.*cos(the0);
                % % vv=abs(vv);
                % [y1,x1]=hist(vv,50);
                % figure;
                % y1=y1/sum(y1);
                % bar(x1,y1);
                % xlim([-10,10]);
                % sum(y1.*x1)
                
                % figure %plot the motion of particles (part 1)
                % h(1)=axes('position',[0.1 0.1,0.8,0.8]);
                % contour(h(1),x,y,phi)
                % h(2)=axes('position',[0.1 0.1,0.8,0.8]);
                % set(h(2),'color','none')
                % linkaxes(h,'xy')
                % % %save video
                % % filename='test_shipin1';
                % % writerObj = VideoWriter(filename);
                % % writerObj.FrameRate = 30;
                % % open(writerObj);
                
                dlx_net=-3*dl:dl/10:3*dl;
                dly_net=-3*dl:dl/10:3*dl;
                [dlx_net,dly_net]=meshgrid(dlx_net,dly_net);
                fxy=1/4/dl^2*exp(-pi/4/dl^2*(dlx_net.^2+dly_net.^2));
                dlxy(:,:,1)=dlx_net;
                dlxy(:,:,2)=dly_net;
                
                kk=0;
                zz=ones(1,N);%number of particles escaped
                clear N1;
                clear t_count;
                for ii=t(2:end)
                    s=sprintf('Ratio=%.1f, num=%d, Pe=%d, n=%.2f, %.4f',rat(kk_r), nnum, PPe, nn, ii/t(end)+ll);
                    disp(s);
                    kk=kk+1;
                    t_count(kk)=ii;
                    for jj=1:N
                        if ismember(nc(jj),nc_boundry)
                            zz(jj)=0;
                        end
                    end
                    N1(kk)=sum(zz);
                    
                    % %plot the motion of particles (part 2)
                    % plot(h(2),x0,y0,'o','markersize',5)
                    % set(h(2),'color','none')
                    % xlim([x(1),x(end)]);
                    % ylim([y(1),y(end)]);
                    % title(['t=',num2str(ii)])
                    % pause(0.1)
                    % % frame =  getframe(gcf);
                    % % writeVideo(writerObj,frame);
                    
                    %random velocity
                    dxyr=dist(2,N,dlxy,fxy);
                    dxr=dxyr(1,:);
                    dyr=dxyr(2,:);
                    
                    %motion + collision
                    x1=x0;
                    y1=y0;
                    index_0=find(nc~=0);%in the case of inside cells
                    ww=vv(index_0)./r0(index_0);
                    dthe=ww*dt;
                    the1=the0(index_0)+dthe;
                    x1(index_0)=x_B0(nc(index_0))+r0(index_0).*cos(the1)+dxr(index_0);
                    y1(index_0)=y_B0(nc(index_0))+r0(index_0).*sin(the1)+dyr(index_0);
                    
                    index_1=find(nc==0);%in the case of out of cells
                    vx=vv(index_1).*sin(the0(index_1));
                    vy=-vv(index_1).*cos(the0(index_1));
                    x1(index_1)=x0(index_1)+dxr(index_1)+vx*dt;
                    y1(index_1)=y0(index_1)+dyr(index_1)+vy*dt;
                    
                    x1(x1<min(x_net(1,:)))=min(x_net(1,:));
                    x1(x1>max(x_net(1,:)))=max(x_net(1,:));
                    y1(y1<min(y_net(:,1)))=min(y_net(:,1));
                    y1(y1>max(y_net(:,1)))=max(y_net(:,1));
                    
                    % reset cells, radii, angles
                    [nc,r0,the1]=find_circle(N,x1,y1,x0,y0,x_B0,y_B0,RR0,the0,nc,vv,dt);
                    the0=the1;
                    x0=x1;
                    y0=y1;
                end
                
                logN=log(N1/N);
                index=find(t_count>=2*tau_eff);
                p=polyfit(t_count(index),logN(index),1);
                logN2=polyval(p,t_count(index));
                tau=-1/p(1);
                D_eff_tem(ll)=((nnum-1)/2*l)^2/tau;
                
            end
            % figure;
            % plot(t_count,logN,'linewidth',2);
            % hold on
            % plot(t_count(index),logN2,'color',[0,0.5,0],'linewidth',2)
            % xlabel('t');
            % ylabel('ln(N/N0)')
            
        end%finish scanning n
        
        D_eff0(kk_num,kk_r)=sqrt(D0*v_p*d);
        D_eff(kk_num,kk_r)=mean(D_eff_tem);
        ddt(kk_n)=log(v_p*dt/d)/(log(10));
        Pe(kk_n)=PPe;
        R0(kk_num,kk_r)=RR;
        wid0(kk_num,kk_r)=wwid;
        
    end%finish scanning num
    
end%finish scanning radio

% % close(writerObj);
eval(['save(''/Users/zhangby/Desktop/cell transport/result/D_eff_number_',num2str(N),'_m_',num2str(m),'.mat'',''D_eff'',''D_eff0'',''d'',''D0'',''Pe'',''n'',''m'',''N'',''rat'',''num'',''R0'',''wid0'')'])




%% check cells
% figure;
% ax(1)=axes('position',[0.1,0.1,0.8,0.8]);
% quiver(ax(1),x_net,y_net,Bx,By);
% ax(2)=axes('position',[0.1,0.1,0.8,0.8]);
% contour(ax(2),x,y,phi)
% set(ax(2),'color','none')
% linkaxes(ax,'x')
% linkaxes(ax,'y')

%Pe=v0*d/D
