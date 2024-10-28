clc;            
clear all;      %Clears the screen
xmax=35;
tmax=0.2;
nx=100000;
nt=100;
x=linspace(0.1,xmax,nx);
t=linspace(0,tmax,nt);
dx=x(2)-x(1);
dt=t(2)-t(1);
p=zeros(nt,nx); %Phi
P=zeros(nt,nx); %Pi
a=zeros(nt,nx); %gxx
M=zeros(nt,nx); %mass
N=zeros(nt,nx); %alpha (lapse)

dpdx=zeros(nt,nx); %Phi'
dadx=zeros(nt,nx); %a'
dNdx=zeros(nt,nx); %alpha'

%Holonomy parameter`
k=0;

%Initial conditions
phi0=10^(-6);
x0=25;
d=1.5;
q=2;
w=1;
for i=1:nx
    p(1,i)=phi0*(w*x(i)^(w-1)-x(i)^w*q*((x(i)-x0)/d)^(q-1)/d)*exp(-((x(i)-x0)/d)^q);%phi=phi0*x^w*exp(-((x-x0)./d).^q)
    P(1,i)=0;
end
for l=1:nx-1
    i=nx+1-l;
    dpdx(1,i)=(p(1,i)-p(1,i-1))/dx;%initial p'
end
%Boundary conditions
for j=1:nt
    a(j,nx)=1;
    M(j,nx)=(1-1/a(j,nx)^2)*x(nx)/2;%Mass
    dadx(j,nx)=0;%Boundary a'=0
    N(j,nx)=1;%Boundary N'=0 already taken
end

%%% Finite difference method %%%
% Time cycle %
for j=1:nt-1
    % Radial cycle %
    for l=1:nx-1
        i=nx+1-l;

        % k1
        k1_dadx = -a(j,i)*((a(j,i)^2-1)/(2*x(i)) - (x(i)/2)*(P(j,i)^2*(1+k^2/a(j,i)^2) + p(j,i))^2);
        
        % Intermediate values
        a_mid = a(j,i) - (dx/2) * k1_dadx;
        k2_dadx = -a_mid*((a_mid^2-1)/(2*x(i)) - (x(i)/2)*(P(j,i)^2*(1+k^2/a_mid^2) + p(j,i))^2);
        
        % Update 
        a(j,i-1) = a(j,i) - dx * k2_dadx;
        dadx(j,i-1) = k2_dadx;

        % RK2 for N equation
        % k1
        k1_dNdx = -N(j,i)*((1 - a(j,i)^2)/(2*x(i)) - (x(i)/2)*(P(j,i)^2*(1 + 3*k^2/a(j,i)^2) + p(j,i))^2);

        % Intermediate values
        N_mid = N(j,i) - (dx/2) * k1_dNdx;
        k2_dNdx = -N_mid*((1 - a_mid^2)/(2*x(i)) - (x(i)/2)*(P(j,i)^2*(1 + 3*k^2/a_mid^2) + p(j,i))^2);

        % Update N 
        N(j,i-1) = N(j,i) - dx * k2_dNdx;
        dNdx(j,i-1) = k2_dNdx;

        % Update M using new value of a
        M(j,i-1) = (1 - 1/a(j,i-1)^2) * x(i-1) / 2;
    end
    
    % P,p EoMs remain as is
    for l=1:nx-1
        i=nx+1-l;
        p(j+1,i) = p(j,i) + (dt/dx)*(N(j,i)/a(j,i)*P(j,i)*(1 + k^2/a(j,i)^2) - N(j,i-1)/a(j,i-1)*P(j,i-1)*(1 + k^2/a(j,i-1)^2));
        P(j+1,i) = P(j,i) + dt/(x(i)^2*dx)*(x(i)^2*N(j,i)/a(j,i)*p(j,i) - x(i-1)^2*N(j,i-1)/a(j,i-1)*p(j,i-1));
        dpdx(j+1,i-1) = (p(j+1,i) - p(j+1,i-1)) / dx;
    end
end


N(N>5)=nan;
N(N<0)=nan;

[T,R]=meshgrid(t,x);

figure(1)
mesh(x,t,N)
%zlim(ax1,[0 8])
title('lapse (N)');
xlabel('x');
ylabel('t');
zlabel('N');

M(M>2)=nan;
M(M<-0.5)=nan;


figure(2)
mesh(x,t,M);
title('mass (M)');
xlabel('x');
ylabel('t');
zlabel('M');


p(p>10)=nan;
p(p<-10)=nan;

figure(3)
mesh(x,t,p);
title('Scalar field \Phi');
xlabel('x');
ylabel('t');
zlabel('\Phi');


P(p>10)=nan;
P(p<-10)=nan;

figure(4)
mesh(x,t,P);
title('Scalar field P_\phi');
xlabel('r');
ylabel('t');
zlabel('\Pi');