
clear all

r=10;
A=pi*r^2;

a=linspace(r,3*r,100);
b=A/pi*ones(size(a))./a;

%%
for ii=1:numel(a)
    p(ii)=ellipsePerimeter([a(ii) b(ii)]);
    ecc(ii)=axes2ecc([a(ii) b(ii)]);
    rain(ii,:)=rainbow(1,1,1,ii,numel(a));
     [la{ii},lo{ii}] = ellipse1(0,0,[a(ii),ecc(ii)]);   
     d(ii)=max(la{ii})-min(la{ii});
end
%%
iq=4*pi*A./p.^2;
iqs=iq-min(iq);
iqs=round(iqs/max(iqs)*(numel(a)-1))+1;

figure(1000)

subplot(211)
 hold on
for ii=1:numel(a)   
   
    plot(la{ii},lo{ii},'color',rain(iqs(ii),:));
    axis equal tight  
end
colormap(rain);
cb=colorbar;
yl=get(cb,'ylim');
set(cb,'ylim',[1 100]);
yt=[1:9:100];
set(cb,'ytick',yt);
ytl=cellfun(@(x) sprintf('%2.2f',x),  num2cell(iq(fliplr(yt))),'uniformoutput',false);
set(cb,'yticklabel',ytl);

%%


ch=4/pi*A*ones(size(a))./d.^2
iqs=ch-min(ch);
iqs=round(iqs/max(iqs)*(numel(a)-1))+1;

figure(1000)
subplot(212)
 hold on
for ii=1:numel(a)      
    plot(la{ii},lo{ii},'color',rain(iqs(ii),:));
    axis equal tight  
end
colormap(rain);
cb=colorbar;
yl=get(cb,'ylim');
set(cb,'ylim',[1 100]);
yt=[1:9:100];
set(cb,'ytick',yt);
ytl=cellfun(@(x) sprintf('%2.2f',x),  num2cell(ch(fliplr(yt))),'uniformoutput',false);
set(cb,'yticklabel',ytl);





iq=4*pi*A./p.^2;


%%

figure(1)
plot(a)
hold on
plot(b,'r')
figure(10)
plot(a/r,p,'black')
figure(100)
plot(a/r,iq,'black')
hold on
plot(a/r,sqrt(iq),'r')

%%
figure(1e5)
clf
N=1000;
M=3;
sss=5;
n=8;
NRA=ones(1,M);
RA =logspace(0,1,M);
for ii=1:M-2
    nonrandbit=NRA(ii);
    randbit=RA(ii);
    t=linspace(0,2*pi,n);
    tt=linspace(0,2*pi,N);
    r=nonrandbit+randbit*rand(size(t));
    %r=nonrandbit+randbit*ones(size(t));
    r(end)=r(1);
    rr=spline(t,[0 r 0],tt);
    x=rr.*cos(tt);
    y=rr.*sin(tt);
    % then use arclength parameterization
    s = zeros(N,1);
    s(1) = 0;
    s(2:N) = sqrt(diff(x).^2+diff(y).^2);
    s = cumsum(s);
    % L == s(N) == sum(sqrt(diff(x).^2+diff(y).^2));
    % interpolation again
    cs = spline(s, [x;y]);
    ss = linspace(s(1),s(N),N);
    ds = (s(N)-s(1))/N;
    curve = ppval(cs,ss);
    curve = curve';
    xx = curve(:,1);
    yy = curve(:,2);
    Area=polyarea(xx,yy);
    P=sum(abs(hypot(diff(xx),diff(yy))));
    [A,B]=meshgrid(xx,yy);
    tmp=hypot(A-A',B-B');
    d=max(tmp(:));    
    iq=4*pi*Area./P.^2;
    ch=4/pi*Area./d.^2;
    leg{ii}=sprintf('iq: %1.2f, ch2: %1.2f',iq,ch);
    hold on
    strtch=1+ii/M*2;
    plot(xx,yy,'color',rainbow(0,1,1,ii,M),'linewidth',2)
end


for ii=M-1:M
    nonrandbit=NRA(ii);
    randbit=RA(ii);
    t=linspace(0,2*pi,n);
    tt=linspace(0,2*pi,N);
    r=nonrandbit+randbit*rand(size(t));
    %r=nonrandbit+randbit*ones(size(t));
    r(end)=r(1);
    rr=spline(t,[0 r 0],tt);
    x=rr.*cos(tt);
    y=1/sss*rr.*sin(tt);
    % then use arclength parameterization
    s = zeros(N,1);
    s(1) = 0;
    s(2:N) = sqrt(diff(x).^2+diff(y).^2);
    s = cumsum(s);
    % L == s(N) == sum(sqrt(diff(x).^2+diff(y).^2));
    % interpolation again
    cs = spline(s, [x;y]);
    ss = linspace(s(1),s(N),N);
    ds = (s(N)-s(1))/N;
    curve = ppval(cs,ss);
    curve = curve';
    xx = curve(:,1);
    yy = curve(:,2);
    Area=polyarea(xx,yy);
    P=sum(abs(hypot(diff(xx),diff(yy))));
    [A,B]=meshgrid(xx,yy);
    tmp=hypot(A-A',B-B');
    d=max(tmp(:));
    iq=4*pi*Area./P.^2;
    ch=4/pi*Area./d.^2;
    leg{ii}=sprintf('iq: %1.2f, ch2: %1.2f',iq,ch);
    hold on
    
    

    plot(xx,yy,'color',rainbow(0,1,1,ii,M),'linewidth',2)

    
%     strtch=1+ii/M*2;
%     xxo=xx-min(xx);
%     xxo=xxo/max(xxo)*strtch;
%     yyo=yy-min(yy);
%     yyo=yyo/max(yyo)*strtch;
%     plot(xxo-mean(xxo),yyo-mean(yyo),'color',rainbow(1,1,1,ii,M))
end
axis equal off tight
legend(leg)
set(gca,'xtick',[])
set(gca,'ytick',[])
%   savefig('../PLOTS/',100,1200,800,['iq2ch2B'])
% savefig('../PLOTS/',100,800,600,['iq2ch2A'])
 