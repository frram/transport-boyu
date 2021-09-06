function [dxy]=dist(varargin)
%scale2: f=nf=n=x00;
%scale1: Nsx=x,Nfx=x0=x-1
%输入：dist(dim,N,x0,f)
%返回按照distribution所取的随机数



if ~isempty(varargin)
    dim=varargin{1};
    if dim==1
        N0=1:varargin{2};%点数
        x0=varargin{3};%网格中心
        dx0=mean(diff(x0));
        x=[x0(1)-dx0/2,x0+dx0/2];%网格边界
        f=varargin{4};%分布
    elseif dim==2
        N0=1:varargin{2};
        xy00=varargin{3};
        x00=xy00(:,:,1);
        y00=xy00(:,:,2);
        x0=x00(1,:);
        y0=y00(:,1)';
        dx0=mean(diff(x0));
        x=[x0(1)-dx0/2,x0+dx0/2];
        dy0=mean(diff(y0));
        y=[y0(1)-dy0/2,y0+dy0/2];
        f=varargin{4};
    end
else
    dim=1;
    N0=1:500;
    x=-3:0.3:3;
    x0=(x(1:end-1)+x(2:end))/2;
    f=exp(-(x0).^2/2);
%     f=sin(x0)+1;
    dx0=mean(diff(x0));
    
%     dim=2;
%     N0=1:10000;
%     x=-2:0.04:2;
%     y=-4:0.04:4;
%     x0=(x(1:end-1)+x(2:end))/2;
%     y0=(y(1:end-1)+y(2:end))/2;
%     x00=meshgrid(x0,y0);
%     y00=meshgrid(y0,x0)';
%     f=exp(-((x00).^2+(y00/2).^2)/2);
%     dx0=mean(diff(x0));
%     dy0=mean(diff(y0)); 
end

if dim==1
    nf=round(f/min(f));
    Nf=sum(nf);
    Ns=zeros(size(x));
    Ns(1)=0;
    for ii=2:length(x)
        Ns(ii)=Ns(ii-1)+nf(ii-1);
    end
    n=zeros(size(x0));
    dx=zeros(size(N0));
    for ii=N0
        r=rand(1);
        R=round(r*Nf);
        index=find(Ns>=R,1);
        if index==1
            index=2;
        end
        n(index-1)=n(index-1)+1;
        dx(ii)=x(index-1)+r*dx0;
    end
    dxy=dx;
elseif dim==2
    nf=round(f/min(min(f)));
    Nfx=sum(nf);
    Nsx=zeros(size(x));
    Nsx(1)=0;
    for ii=2:length(x)
        Nsx(ii)=Nsx(ii-1)+Nfx(ii-1);
    end
    Nsy=zeros([size(y00,1)+1,size(y00,2)]);
    Nfy=sum(nf);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nsy(1,:)=0;
    for ii=1:length(x0)
        for jj=2:length(y)
            Nsy(jj,ii)=Nsy(jj-1,ii)+nf(jj-1,ii);
        end
    end
    n=zeros(size(x00));
    dx=zeros(size(N0));
    dy=zeros(size(N0));
    for ii=N0
        rx=rand(1);
        Rx=round(rx*sum(Nfx));
        index_x=find(Nsx>=Rx,1);
        if index_x==1
            index_x=2;
        end
        
        ry=rand(1);
        Ry=round(ry*Nfy(index_x-1));
        index_y=find(Nsy(:,index_x-1)>=Ry,1);
        if index_y==1
            index_y=2;
        end
        n(index_y-1,index_x-1)=n(index_y-1,index_x-1)+1;
        dx(ii)=x(index_x-1)+rx*dx0;
        dy(ii)=y(index_y-1)+ry*dy0;
    end
    dxy=[dx;dy];
end




% %检验
% if dim==1
% %     figure
% %     subplot(2,1,1)
% %     plot(x0,f);
% %     subplot(2,1,2)
% %     bar(x0,n)
% %     %%%检验
% %     [y1,x1]=hist(dx,100);
% %     y1=y1/length(dx)/mean(diff(x1));
% %     figure;
% %     bar(x1,y1);
% 
% elseif dim==2
%     figure;
%     subplot(2,1,1)
%     surf(x00,y00,f);
%     shading interp
%     subplot(2,1,2)
%     surf(x00,y00,n)
%     shading interp
% end
% 1