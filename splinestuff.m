e=load('e.mat')
%%
dx=10000;
ssh=e.AntiCycs(1).profiles.x.ssh
 x=linspace(-1,1,numel(ssh))*dx
 X=linspace(-1,1,numel(ssh)*10)*dx
  SSH=spline(x,ssh,X)
  %%
  
Fs = 1/dx;                    % Sampling frequency
T = 1/Fs;                     % Sample dist
L = numel(X);                 % ?
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(SSH,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
ff=[-fliplr(f) -fliplr(f)]
%%
% Plot single-sided amplitude spectrum.
clf
loglog(f,2*abs(Y(1:NFFT/2+1))) 
xlabel('Frequency (1/m)')
ylabel('|Y(f)|')
xt=get(gca,'xtick')
set(gca,'xticklabel',xt.^(-1))
%%
clf
plot(abs(Y))
xlabel('Frequency (1/m)')
ylabel('|Y(f)|')
xt=get(gca,'xtick')
set(gca,'xticklabel',xt.^(-1))
%%
Y=fft(SSH)
Y(1)=0
Y(end)=0
Y(20:end-20)=0
Y=smooth(Y,10)
iY=ifft(Y)
subplot(311)
plot(abs(Y))
subplot(312)
plot(abs(iY))
subplot(313)
plot(SSH)
%%
figure(2)
subplot(221)
plot(X,SSH)
subplot(222)
SSHd=detrend(SSH)
plot(X,SSHd)
dS=diff(SSHd,1)*1e6
dS=spline((X(2:end)+X(1:end-1))/2,dS,X)
subplot(223)
plot(dS)
grid on
dSS=diff(smooth(SSHd,3),2)
dSS=[nan; dSS; nan]
subplot(224)
plot(dSS)
grid on
%%
dSd=diff(sign(dS))
xx=find(dSd' & dSS(1:end-1)>0)
xl.a=xx(1);
xl.b=xx(2)+1;
%%
clf
ssho=ssh;
xo=x;
ssh=(SSH);
x=X;

[f,gof,out] = fit(x.',detrend(ssh)','fourier2')
%%
[f2]=diffCentered(2,x,f(x));
[f2r]=diffCentered(2,x,ssh);
[f1]=diffCentered(1,x,f(x));
[f1r]=diffCentered(1,x,ssh);
clf
hold on
a(1)=plot(x,(ssh'),'b--')
a(2)=plot(x,(f(x)),'b')
%plot(x,normc(f1'),'black')
a(3)=plot(x,normc(f2r'),'r--')
a(4)=plot(x,normc(f2'),'r')
plot(xo,(ssho'),'*','markersize',8)
%plot(x,normc(f1r'),'black--')
grid on
legend('ssh','first diff','second diff')
legend('intrp(ssh)','Fourier2(detrend(#1))','second diffs([#1 #2])')
set(gca,'ytick',0,'xtick',[])
axis([-1e4 1e4 -.6 .4])
axis tight
set(a(:),'linewidth',2)
set(get(gcf,'children'),'linewidth',2)
xlabel('all normalized')


%%
clf
plot(normc(diff(f(X))))
hold on
plot(normc(diff(f(X),2)),'r')
hold on
plot(normc(diff(SSH,2)'),'r--')
% hold on
% plot(normc(SSHd'),'g')
hold on
plot(normc(SSH'),'y')
hold on
plot(normc(f(X)),'b--')
grid on





