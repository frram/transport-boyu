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
N=2000;
d=1;
D0=1;
tau0=d^2/D0;

% Pe=[25,50,100,400,900];
Pe = [25];
n=[1,2,5,10,20,50,80,100];

%for saving data
data.D_eff1=zeros(length(Pe)+1,length(n)+1); %simulation results with method1
data.D_eff1(1,2:end)=n;
data.D_eff1(2:end,1)=Pe;
data.D_eff2=zeros(length(Pe)+1,length(n)+1); %simulation results with method2
data.D_eff2(1,2:end)=n;
data.D_eff2(2:end,1)=Pe;
data.D_eff0=zeros(length(Pe)+1,length(n)+1);%theoretical results
data.D_eff0(1,2:end)=n;
data.D_eff0(2:end,1)=Pe;

% set magnetic cells: centers and radii
x_B0=-9.5*d:d:10.5*d;
y_B0=0.5*d;
[x_B0,y_B0]=meshgrid(x_B0,y_B0);
x_B0=x_B0(:)';
y_B0=y_B0(:)';
n0=length(x_B0);
R0=ones(1,n0)*d/2.01;

%set phi and magnetic field #set like Rosenbluth
x=x_B0(1)-0.5*d:0.01*d:x_B0(end)+0.5*d;
y=y_B0(1)-0.5*d:0.01*d:y_B0(end)+0.5*d;
[x,y]=meshgrid(x,y);
phi=sin(pi*x/d).*sin(pi*y/d);
[By,Bx]=grad(phi,x,y);
x_net=x(2:end-1,2:end-1);
y_net=y(2:end-1,2:end-1);
By=By;
Bx=-Bx;

% %plot the magnetic field
% figure
% set(gcf,'position',[100,50,600,550])
% h(1)=axes('position',[0.1 0.1,0.8,0.8]);
% contour(h(1),x,y,phi)
% hold on
% plot(h(1),x_B0,y_B0,'o','color','r')
% h(2)=axes('position',[0.1 0.1,0.8,0.8]);
% quiver(h(2),x_net,y_net,Bx,By)
% set(h(2),'color','none')
% xlim([x_B0(1)-0.5*d,x_B0(end)+0.5*d])
% ylim([y_B0(1)-0.5*d,y_B0(end)+0.5*d])
% linkaxes(h,'xy')

% scan Pe
data1=[];
data2=[];
kk_p=0;
for PPe=Pe
    kk_p=kk_p+1;
    v_p=PPe*D0/d;
    
    %scan dt
    kk_n=0;
    for nn=n
        kk_n=kk_n+1;
        
        dt=d/v_p/nn;
        tau_eff=tau0/sqrt(PPe);
        t=0:dt:5*tau_eff;
        sig=d/sqrt(PPe*nn);
        mu=0;
        sigma=sig^2;
        
        for ll=1:2%averaging
            %initialize the particle coordinates
            nc0=(n0+1)/2*ones(1,N);
            nc=nc0;
            x0=x_B0(nc)+rand(1,N)*d-0.5*d;
            y0=y_B0(nc)+rand(1,N)*d-0.5*d;
            
%             figure %check the initial particles
%             h(1)=axes('position',[0.1 0.1,0.8,0.8]);
%             contour(h(1),x,y,phi)
%             hold on
%             plot(h(1),x_B0,y_B0,'o','color','r')
%             h(2)=axes('position',[0.1 0.1,0.8,0.8]);
%             plot(h(2),x0,y0,'o')
%             set(h(2),'color','none')
%             linkaxes(h,'xy')
            
            posi=zeros(N,length(t),2);%position of every particle
            posi(:,1,1)=x0;
            posi(:,1,2)=y0;
            posi_nc=ones(N,length(t))*(n0+1)/2;%cells that particles locate in
%             vv=v_p*ones(1,N);%1，固定速度
            vv=mvnrnd(0,(v_p*sqrt(2*pi)/2)^2,N)';%set parallel velocity distribution
            
            % figure%plot the motion of particles (part 1)
            % h(1)=axes('position',[0.1 0.1,0.8,0.8]);
            % contour(h(1),x,y,phi,20)
            % h(2)=axes('position',[0.1 0.1,0.8,0.8]);
            % set(h(2),'color','none')
            % linkaxes(h,'xy')
            
            kk=0;
            for ii=t(2:end)
                s=sprintf('%d,%.2f,%.4f',PPe,nn,ii/t(end)+ll);
                disp(s);
                kk=kk+1;
                t_count(kk)=ii;
                
                % %plot the motion of particles (part 2)
                % plot(h(2),x0,y0,'o','markersize',5)
                % set(h(2),'color','none')
                % xlim([-5*d,5*d]);
                % ylim([y(1),y(end)]);
                % pause(0.1)
                
                % random velocity
                dxyr=mvnrnd(mu,sigma,N);
                dxr=dxyr(:,1);
                
                %motion + collision: Rosenbluth cell
                r0 = [x0,y0]'; 
                tspan = [0:dt/2:dt];
                f = @(t,y) [-vv'.*sin(y(1:N)*pi/d).*cos(y(N+1:2*N)*pi/d); vv'.*cos(y(1:N)*pi/d).*sin(y(N+1:2*N)*pi/d)];
                [tt,r] = ode45(f,tspan,r0,odeset('RelTol',1e-7));
                x0=r(3,1:N)+dxr';
                y0=r(3,N+1:2*N);
                
                nc=floor(x0)+nc0;
                posi(:,kk+1,1)=x0;
                posi(:,kk+1,2)=y0;
                posi_nc(:,kk+1)=nc;
                
                % % motion + collision: round cell
                % x1=x0;
                % y1=y0;
                % index_0=find(nc~=0);%inside cells
                % ww=vv(index_0)./r0(index_0);
                % dthe=ww*dt;
                % the1=the0(index_0)+dthe;
                % x1(index_0)=x_B0(nc(index_0))+r0(index_0).*cos(the1)+dxr(index_0)';
                % y1(index_0)=y_B0(nc(index_0))+r0(index_0).*sin(the1);
                %
                % index_1=find(nc==0);%outside of cells
                % vx=vv(index_1).*sin(the0(index_1));
                % vy=-vv(index_1).*cos(the0(index_1));
                % x1(index_1)=x0(index_1)+vx*dt+dxr(index_1)';
                % y1(index_1)=y0(index_1);
                %
                % x1(x1<min(x_net(1,:)))=min(x_net(1,:));
                % x1(x1>max(x_net(1,:)))=max(x_net(1,:));
                % y1(y1<min(y_net(:,1)))=min(y_net(:,1));
                % y1(y1>max(y_net(:,1)))=max(y_net(:,1));
                % %reset cells, radii, angles
                % [nc,r0,the0]=find_circle(N,x1,y1,x0,y0,x_B0,y_B0,R0,the0,nc,vv,dt);
                % x0=x1;
                % y0=y1;
            end
            
            n_interval=floor((length(t)-1)/2);
            t_interval=dt:dt:dt*n_interval;
            
            %Method 1：calculate MSD
            MSD_tpt=zeros(1,n_interval,2);
            for ii = 1:n_interval
                s1=sprintf('method1,%.3f',ii/n_interval+ll);
                disp(s1)
                MSD_tpt(:,ii,1)=mean( mean(  (  posi(:,(1:n_interval)+ii,1) - posi(:,(1:n_interval),1)  ).^2, 2    ));
                MSD_tpt(:,ii,2)=mean( mean(  (  posi(:,(1:n_interval)+ii,2) - posi(:,(1:n_interval),2)  ).^2, 2    ));
            end
            MSD=sum(MSD_tpt,3);
            MSD=MSD_tpt(:,:,1);%1-D case
            
            index=find(t_interval>=tau_eff);
            % index=find(t_interval>=0);
            t_interval1=t_interval(index);
            MSD1=MSD(index);
            p=polyfit(t_interval1,MSD1,1);
            MSD_fit=polyval(p,t_interval1);
            D1=p(1);
            data_tpt=[PPe;nn;D1];
            data1=[data1,data_tpt]
%             %plot to check linear or not
%             figure;
%             plot(t_interval,MSD,'bo','markersize',5);
%             hold on
%             plot(t_interval1,MSD_fit,'r','linewidth',2);
%             xlabel('\Deltat')
%             ylabel('MSD')
%             legend({'simulation','fitting'},'fontsize',20)
%             title(['dt=',num2str(dt),', d=',num2str(d)])
            
            
            %Method 2：count how many cells that particles have transported
            MSD=zeros(1,n_interval);
            for ii = 1:n_interval%different t_interval
                s2=sprintf('method1,%.3f',ii/n_interval+ll);
                disp(s2)
                MSD_tpt=0;
                for jj = 1:n_interval%averaging for different t_start
                    MSD_tpt = MSD_tpt+sum(((posi_nc(:,jj)-posi_nc(:,jj+ii))*d).^2);
                end
                MSD(ii)=MSD_tpt/N/n_interval;
            end
            index=find(t_interval>=tau_eff);
            t_interval1=t_interval(index);
            MSD2=MSD(index);
            p=polyfit(t_interval1,MSD2,1);
            MSD_fit=polyval(p,t_interval1);
            D2=p(1);
            data_tpt=[PPe;nn;D2];
            data2=[data2,data_tpt]
%             %plot to check linear or not
%             figure;
%             plot(t_interval,MSD,'bo','markersize',5);
%             hold on
%             plot(t_interval1,MSD_fit,'r','linewidth',2);
%             xlabel('\Deltat')
%             ylabel('MSD')
%             legend({'simulation','fitting'},'fontsize',20)
%             title(['dt=',num2str(dt),', d=',num2str(d)])
        end
        data.D_eff1(kk_p+1,kk_n+1)=mean(data1(3,(kk_p-1)*length(n)*ll+(kk_n-1)*ll+1:(kk_p-1)*length(n)*ll+(kk_n-1)*ll+ll));
        data.D_eff2(kk_p+1,kk_n+1)=mean(data2(3,(kk_p-1)*length(n)*ll+(kk_n-1)*ll+1:(kk_p-1)*length(n)*ll+(kk_n-1)*ll+ll));
        data.D_eff0(kk_p+1,kk_n+1)=sqrt(D0*v_p*d);
    end %finish scanning n
    
end %finish scanning Pe

% figure;plot(data.D_eff(1,2:end),data.D_eff(2:end,2:end),'o',data.D_eff0(1,2:end),data.D_eff0(2:end,2:end))
% figure;plot(n,data.D_eff(2:end,2:end),'o',n,data.D_eff0(2:end,2:end))


