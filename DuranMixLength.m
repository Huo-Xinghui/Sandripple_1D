% Smooth surface
ustar=0.9;
nu=1.51e-5;
kapa=0.4;
% c=9.98e-4;
% Rc=8.27;
c=0;
Rc=5;
num=49;
zmax=0.3;
zctrl=1.2;
n=linspace(0,zmax,num);
z=zeros(1,num+1);
dz=zeros(1,num);
z(1)=0;
for i=1:num
    dz(i)=zctrl^(i-1)*zmax/((zctrl^num-1)/(zctrl-1));
    z(i+1)=z(i)+dz(i);
end
zplus=z*ustar/nu;
l1=zeros(1,num+1);
l2=zeros(1,num+1);
dl2=zeros(1,num);
du1=zeros(1,num);
du2=zeros(1,num);
u1=zeros(1,num+1);
u2=zeros(1,num+1);
l1(1)=0;
l2(1)=l1(1);
% l2(1)=(2*kapa*ustar*z(1)^1.5+3*sqrt(Rc)*nu*c)^2/(36*Rc*nu^2);
u1(1)=0;
u2(1)=u1(1);
for i=1:num
    dl2(i)=kapa*(1-exp(-sqrt(1/Rc*u2(i)*l2(i)/nu)));
    l1(i+1)=kapa*z(i+1)*(1-exp(-1/26*zplus(i+1)));
    if i==1
        l2(i+1)=(2*kapa*ustar*z(i+1)^1.5+3*sqrt(Rc)*nu*c)^2/(36*Rc*nu^2);
    else
        l2(i+1)=l2(i)+dl2(i)*dz(i);
    end
    l1c=(l1(i)+l1(i+1))/2;
    l2c=(l2(i)+l2(i+1))/2;
    du1(i)=(-nu+sqrt(nu^2+4*l1c^2*ustar^2))/(2*l1c^2);
    du2(i)=(-nu+sqrt(nu^2+4*l2c^2*ustar^2))/(2*l2c^2);
    u1(i+1)=du1(i)*dz(i)+u1(i);
    u2(i+1)=du2(i)*dz(i)+u2(i);
end
u1plus=u1/ustar;
u2plus=u2/ustar;
% plot(zplus,l1,'k')
% hold on
% plot(zplus,l2,'r')
% hold off
figure(1);
loglog(zplus,u1plus,'k--o')
hold on
loglog(zplus,u2plus,'r')
hold off
savefig(1,"loglog")
figure(2);
semilogx(zplus,u1plus,'k--o')
hold on
semilogx(zplus,u2plus,'r')
hold off
savefig(1,"semilog")