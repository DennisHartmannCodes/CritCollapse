clc;            
clear all;      %Clears the screen
rmax=10;
tmax=5;
nr=1000;
nt=100;
r=linspace(0.1,rmax,nr);
t=linspace(0,tmax,nt);
dr=r(2)-r(1);
dt=t(2)-t(1);
p=zeros(nt,nr); %phi
P=zeros(nt,nr); %P_phi
E=zeros(nt,nr); %E^varphi
a=zeros(nt,nr); %alpha (lapse)
%phi(t,r)

dpdr=zeros(nt,nr); %phi'
ddpddr=zeros(nt,nr); %phi''
dEdr=zeros(nt,nr); %(E^varphi)'
dadr=zeros(nt,nr); %alpha'

%Holonomy parameter
k=0;

%Initial conditions
phi0=0;
r0=8;
d=0.1;
q=2;
for i=1:nr
    p(1,i)=-phi0.*r(i).*exp(-((r(i)-r0)./d).^q);
    P(1,i)=0;
end
for l=1:nr-1
    i=nr+1-l;
    dpdr(1,i)=(p(1,i)-p(1,i-1))/dr;%initial phi'
    ddpddr(1,i)=(dpdr(1,i)-dpdr(1,i-1))/dr;%initial phi''
end
%Boundary conditions
for j=1:nt
    E(j,nr)=r(nr);
    dEdr(j,nr)=1;%Boundary E'=1
    a(j,nr)=1;%Boundary a'=0 already taken
end

%%% Finite difference method %%%
% Time cycle %
for j=1:nt-1
    % Radial cycle %
    for l=1:nr-1
        i=nr+1-l;
        % E,a EoMs %
        E(j,i-1)=E(j,i)-dr*dEdr(j,i);
        dEdr(j,i-1)=-E(j,i)*(E(j,i)^2/(2*r(i)^3)-3/(2*r(i))-2*pi*r(i)*(P(j,i)^2/r(i)^4+dpdr(j,i)^2*cos(k*p(j,i))^2));
        a(j,i-1)=a(j,i)-dr*dadr(j,i);
        dadr(j,i-1)=-a(i)*(2/r(i)-E(j,i)^2/(r(i)^3)-dEdr(j,i)/E(j,i));
    end
    % P,p EoMs %
    for l=1:nr-1
        i=nr+1-l;
        P(j+1,i)=P(j,i)+dt*r(i)^2/E(j,i)*((3*a(j,i)*E(j,i)-r(i)*a(j,i)*dEdr(j,i)+dadr(j,i)*E(j,i)*r(i))/E(j,i)*dpdr(j,i)*cos(k*p(j,i))^2+r(i)*a(j,i)*ddpddr(j,i)*cos(k*p(j,i))^2-r(i)*a(j,i)*k*dpdr(j,i)^2*cos(k*p(j,i)).*sin(k*p(j,i)));
        p(j+1,i)=p(j,i)+dt*a(j,i)/(E(j,i)*r(i))*P(j,i);
        dpdr(j+1,i-1)=(p(j+1,i)-p(j+1,i-1))/dr;
        ddpddr(j+1,i-1)=(dpdr(j+1,i)-dpdr(j+1,i-1))/dr;
    end
end

a(a>10)=nan;
a(a<-1)=nan;

[T,R]=meshgrid(t,r);

figure(1)
mesh(r,t,a)
%zlim(ax1,[0 8])
title('lapse (a)');
xlabel('r');
ylabel('t');
zlabel('a');

E(E>400)=nan;
E(E<-1)=nan;

figure(2)
mesh(r,t,E);
title('triad (E^\varphi)');
xlabel('r');
ylabel('t');
zlabel('E^\varphi/r');


p(p>1)=nan;
p(p<-1)=nan;

figure(3)
mesh(r,t,p);
title('Scalar field \phi');
xlabel('r');
ylabel('t');
zlabel('\phi');


P(p>1)=nan;
P(p<-1)=nan;

figure(4)
mesh(r,t,P);
title('Scalar field P_\phi');
xlabel('r');
ylabel('t');
zlabel('\phi');