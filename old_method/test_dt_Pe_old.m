%d: distance between two cell centers
%D0: pure diffusion coefficient
%v_p: v_parallel
%v_r: v_random
%n: n=d/(v_p*dt) indicates step length of parallel velocity
%Pe: Peclet number

%nc: the number of cells in which the particles are located. if nc=0, then
%out of cell.
%r0 & the0: polar coordinates of particles in the cell
%x0 & y0: vertical coordinates of partilces for ploting figure

%% fix d,D0,Pe, scan n (i.e., fix v_p, scan dt)

clear;clc
d=1;
D0=0.1;
tau0=d^2/D0;

Pe=[20,40,60,80,100,150,200,300,400,600,900];
Pe=fliplr(Pe);
n=[1:0.1:3, 3.2:0.2:5, 6:1:16, 18:2:50, 55:5:100, 110:10:150, 165:15:250];
n=n(end:-1:1);

lg10=-log10(n);

%for saving data
data.D_eff=zeros(length(Pe)+1,length(lg10)+1); %simulation results
data.D_eff(1,2:end)=lg10;
data.D_eff(2:end,1)=Pe;
data.D_eff0=zeros(length(Pe)+1,length(lg10)+1);%theoretical results
data.D_eff0(1,2:end)=lg10;
data.D_eff0(2:end,1)=Pe;

% set magnetic cells: centers and radii
x_B0=-2*d:d:2*d;
y_B0=-2*d:d:2*d;
[x_B0,y_B0]=meshgrid(x_B0,y_B0);
x_B0=x_B0(:)';
y_B0=y_B0(:)';
n0=length(x_B0);
R0=ones(1,n0)*d/2.01;

%set phi and magnetic field #round cell
x=-2.5*d:0.01*d:2.5*d;
y=-2.5*d:0.01*d:2.5*d;
[x,y]=meshgrid(x,y);
phi=zeros(size(x));
for ii=1:n0
    index=find((x-x_B0(ii)).^2+(y-y_B0(ii)).^2<=R0(ii)^2);
    phi(index)=cos(pi*(x(index)-x_B0(ii))/d).*cos(pi*(y(index)-y_B0(ii))/d);
    %     phi(index)=cos(pi*x(index)/d).*cos(pi*y(index)/d);
end
[By,Bx]=grad(phi,x,y);
x_net=x(2:end-1,2:end-1);
y_net=y(2:end-1,2:end-1);
By=By;
Bx=-Bx;

%scan Pe
kk_p=0;
for PPe=Pe
    kk_p=kk_p+1;
    v_p=PPe*D0/d;
    
    %scan n
    kk_n=0;
    for nn=n
        dt=d/v_p/nn;
        kk_n=kk_n+1;
        clear t_count;
        clear N1;
        
        v_r=sqrt(D0/dt);
        dl=v_r*dt;
        
        tau_eff=tau0/sqrt(PPe);
        t=0:dt:6*tau_eff;
        
        for ll=1:9%averaging with 9 times
            %initialize the particle coordinates
            N=2000;
            nc=13*ones(1,N);
            r0=rand(1,N).*R0(nc);
            the0=rand(1,N)*2*pi-pi;
            x0=r0.*cos(the0)+x_B0(nc);
            y0=r0.*sin(the0)+y_B0(nc);
            
            % figure %check the initial particles
            % h(1)=axes('position',[0.1 0.1,0.8,0.8]);
            % contour(h(1),x,y,phi)
            % hold on
            % plot(h(1),x_B0,y_B0,'o','color','r')
            % h(2)=axes('position',[0.1 0.1,0.8,0.8]);
            % plot(h(2),x0,y0,'o')
            % set(h(2),'color','none')
            % linkaxes(h,'xy')
            
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
            
            % figure %plot the motion of particles (part 1)
            % h(1)=axes('position',[0.1 0.1,0.8,0.8])
            % contour(h(1),x,y,phi)
            % h(2)=axes('position',[0.1 0.1,0.8,0.8])
            % set(h(2),'color','none')
            % linkaxes(h,'xy')
            
            % set random velocity step length
            dlx_net=-3*dl:dl/10:3*dl;
            dly_net=-3*dl:dl/10:3*dl;
            [dlx_net,dly_net]=meshgrid(dlx_net,dly_net);
            fxy=1/4/dl^2*exp(-pi/4/dl^2*(dlx_net.^2+dly_net.^2));
            dlxy(:,:,1)=dlx_net;
            dlxy(:,:,2)=dly_net;
            
            kk=0;
            zz=ones(1,N); %number of particles escaped
            nc_boundry=[1,2,3,4,5,6,10,11,15,16,20,21,22,23,24,25];
            for ii=t(2:end)
                s=sprintf('%d,%.2f,%.4f',PPe,nn,ii/t(end)+ll);
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
                % pause(0.1)
                
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
                vx=vv(index_1).*sin(the0(index_1));%keep the velocity while in cells
                vy=-vv(index_1).*cos(the0(index_1));
                x1(index_1)=x0(index_1)+dxr(index_1)+vx*dt;%is vx*dt needed here?
                y1(index_1)=y0(index_1)+dyr(index_1)+vy*dt;
                
                x1(x1<min(x_net(1,:)))=min(x_net(1,:));
                x1(x1>max(x_net(1,:)))=max(x_net(1,:));
                y1(y1<min(y_net(:,1)))=min(y_net(:,1));
                y1(y1>max(y_net(:,1)))=max(y_net(:,1));
                % reset cells, radii, angles
                [nc,r0,the1]=find_circle(N,x1,y1,x0,y0,x_B0,y_B0,R0,the0,nc,vv,dt);
                x0=x1;
                y0=y1;
            end
            
            logN=log(N1/N);
            index=find(t_count>=2*tau_eff);
            p=polyfit(t_count(index),logN(index),1);
            logN2=polyval(p,t_count(index));
            tau=-1/p(1);
            D_eff_tem(ll)=(2*d)^2/tau;
            
            % figure;
            % plot(t_count,logN,'linewidth',2);
            % hold on
            % plot(t_count(index),logN2,'color',[0,0.5,0],'linewidth',2)
            % xlabel('t');
            % ylabel('ln(N/N0)')
            
        end % finish 9 calculations
        % ddt(kk_n)=lg10;
        D_eff0(kk_n)=sqrt(D0*v_p*d);
        D_eff(kk_n)=mean(D_eff_tem);
    end %finish scanning n
    % figure;plot(ddt,D_eff,ddt,D_eff0);
    % eval(['save(''C:\Users\tsing\Desktop\french festival\D_eff_',num2str(N),'_nn_',num2str(nn),'_D0_',num2str(D0),'.mat'',''Pe'',''D_eff'',''D_eff0'')'])
    data.D_eff(kk_p+1,2:end)=D_eff;
    data.D_eff0(kk_p+1,2:end)=D_eff0;
end %finish scanning Pe

figure;plot(data.D_eff(1,2:end),data.D_eff(2:end,2:end),'o',data.D_eff0(1,2:end),data.D_eff0(2:end,2:end))
figure;plot(n,data.D_eff(2:end,2:end),'o',n,data.D_eff0(2:end,2:end))


v_p='Pe*D0/d';
dt='d/v_p*10^(lg10)';
nn='10^(-lg10)';
save('D_eff_total_nn.mat','data','d','D0','Pe','lg10','v_p','dt','n','nn')


% check cells
figure;
ax(1)=axes('position',[0.1,0.1,0.8,0.8]);
quiver(ax(1),x_net,y_net,Bx,By);
ax(2)=axes('position',[0.1,0.1,0.8,0.8]);
contour(ax(2),x,y,phi)
set(ax(2),'color','none')
linkaxes(ax,'x')
linkaxes(ax,'y')

%Pe=v0*d/D

