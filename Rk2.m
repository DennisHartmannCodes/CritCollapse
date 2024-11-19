clc;            
clear all;      %Clears the screen
xmax=35;
tmax=0.2;
nx=1000;
nt=1000;
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
phi0=10^(-3);
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
for j = 1:nt - 1
    % Radial cycle %
    for l = 1:nx - 1
        i = nx + 1 - l;

        % RK2 for dadx equation
        k1_dadx = -a(j,i) * ((a(j,i)^2 - 1) / (2 * x(i)) - (x(i) / 2) * (P(j,i)^2 * (1 + k^2 / a(j,i)^2) + p(j,i))^2);
        a_mid = a(j,i) - (dx / 2) * k1_dadx;

        k2_dadx = -a_mid * ((a_mid^2 - 1) / (2 * x(i)) - (x(i) / 2) * (P(j,i)^2 * (1 + k^2 / a_mid^2) + p(j,i))^2);

        a(j,i - 1) = a(j,i) - dx * k2_dadx;

        dadx(j,i - 1) = k2_dadx;

        % RK2 for dNdx equation
        k1_dNdx = -N(j,i) * ((1 - a(j,i)^2) / (2 * x(i)) - (x(i) / 2) * (P(j,i)^2 * (1 + 3 * k^2 / a(j,i)^2) + p(j,i))^2);
        N_mid = N(j,i) - (dx / 2) * k1_dNdx;
        k2_dNdx = -N_mid * ((1 - a_mid^2) / (2 * x(i)) - (x(i) / 2) * (P(j,i)^2 * (1 + 3 * k^2 / a_mid^2) + p(j,i))^2);
        N(j,i - 1) = N(j,i) - dx * k2_dNdx;
        dNdx(j,i - 1) = k2_dNdx;

        % Update M using the new value of a
        M(j,i - 1) = (1 - 1 / a(j,i - 1)^2) * x(i - 1) / 2;
    end
    
    % RK2 for p and P equations
    for l = 1:nx - 2
        i = nx - l;
        
        % k1 for p and P
        k1_p_ahead = (dt/(2 * dx))*(N(j,i + 1)/a(j,i + 1) * P(j,i + 1) * (1 + k^2 / a(j,i + 1)^2) - N(j,i) / a(j,i) * P(j,i) * (1 + k^2 / a(j,i)^2));
        k1_P_ahead = (dt/(2 * x(i)^2 * dx)) * (x(i + 1)^2 * N(j,i + 1) / a(j,i + 1) * p(j,i + 1) - x(i)^2 * N(j,i) / a(j,i) * p(j,i));
        
        k1_p_behind = (dt/(2 * dx))*(N(j,i)/a(j,i) * P(j,i) * (1 + k^2 / a(j,i)^2) - N(j,i - 1) / a(j,i - 1) * P(j,i - 1) * (1 + k^2 / a(j,i - 1)^2));
        k1_P_behind = (dt/(2 * x(i)^2 * dx)) * (x(i)^2 * N(j,i) / a(j,i) * p(j,i) - x(i - 1)^2 * N(j,i - 1) / a(j,i - 1) * p(j,i - 1));

        % Intermediate values for p and P
        p_mid_ahead = p(j,i) + k1_p_ahead / 2;
        P_mid_ahead = P(j,i) + k1_P_ahead / 2;

        p_mid_behind = p(j,i) + k1_p_behind / 2;
        P_mid_behind = P(j,i) + k1_P_behind / 2;
        
        %recalculate p_mid and P_mid between (j, i + 1) (j, i) and (j, i -1) so that it isnt constant

        % k2 for p and P
        k2_p = (dt / (2 * dx)) * (N(j,i + 1) / a(j,i + 1) * P_mid_ahead * (1 + k^2 / a(j,i + 1)^2) - N(j,i - 1) / a(j,i - 1) * P_mid_behind * (1 + k^2 / a(j,i - 1)^2));
        k2_P = (dt / (2 * x(i)^2 * dx)) * (x(i + 1)^2 * N(j,i + 1) / a(j,i + 1) * p_mid_ahead - x(i - 1)^2 * N(j,i - 1) / a(j,i - 1) * p_mid_behind);
        
        % Update p and P with RK2
        p(j + 1, i) = p(j, i) + k2_p;
        P(j + 1, i) = P(j, i) + k2_P;
        
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

figure(4)
mesh(x,t,P);
title('Scalar field P_\phi');
xlabel('r');
ylabel('t');
zlabel('\Pi');

% Define output file for the GIF
outputFile = 'evolution.gif';

% Set up figure
figure;

% Loop through each time slice to create frames for the GIF
for j = 1:nt
    % Plot p vs x at time slice j
    subplot(2, 2, 1);
    plot(x, p(j, :), 'LineWidth', 1.5);
    title(['\Phi (Scalar Field) at t = ' num2str(t(j))]);
    xlabel('x');
    ylabel('\Phi');
    ylim([-10, 10]); % Adjust limits as necessary for your data

    % Plot P vs x at time slice j
    subplot(2, 2, 2);
    plot(x, P(j, :), 'LineWidth', 1.5);
    title(['\Pi (Scalar Field) at t = ' num2str(t(j))]);
    xlabel('x');
    ylabel('\Pi');
    ylim([-10, 10]); % Adjust limits as necessary for your data

    % Plot N vs x at time slice j
    subplot(2, 2, 3);
    plot(x, N(j, :), 'LineWidth', 1.5);
    title(['Lapse (N) at t = ' num2str(t(j))]);
    xlabel('x');
    ylabel('N');
    ylim([0, 5]); % Adjust limits as necessary for your data

    % Plot M vs x at time slice j
    subplot(2, 2, 4);
    plot(x, M(j, :), 'LineWidth', 1.5);
    title(['Mass (M) at t = ' num2str(t(j))]);
    xlabel('x');
    ylabel('M');
    ylim([-0.5, 2]); % Adjust limits as necessary for your data

    % Capture the frame as an image
    frame = getframe(gcf);
    img = frame2im(frame);
    [imgIndexed, cmap] = rgb2ind(img, 256);
    
    % Write to GIF file
    if j == 1
        % For the first frame, create the file
        imwrite(imgIndexed, cmap, outputFile, 'gif', 'LoopCount', Inf, 'DelayTime', .01);
    else
        % For subsequent frames, append to the file
        imwrite(imgIndexed, cmap, outputFile, 'gif', 'WriteMode', 'append', 'DelayTime', .01);
    end
end
