% This program calls the ODE solver to solve the system.
% The spatial discretizaton is determined using the spectral
% method. 

% Initial vector
L=3;
N=32;
sy=8*N;
y = zeros(sy,1);
y(1:N)=0;                        % Vd
y(N+1:2*N)=0;                    % Vsh
y(2*N+1:3*N)=5.01;               % Ca
y(3*N+1:4*N)=750*10^6;           % Rss
y(4*N+1:5*N)=35;                 % nbar

% Initial Vector for the Hodgkin Huxley variables n,m and h

alfn=0.1/(exp(1)-1);
betn=0.125;
alfm=2.5/(exp(2.5)-1);
betm=4;
alfh=0.07;
beth=1/(exp(3)+1);
y(5*N+1:6*N)=alfn/(alfn+betn);   % n
y(6*N+1:7*N)=alfm/(alfm+betm);   % m
y(7*N+1:8*N)=alfh/(alfh+beth);   % h

tic;
[D2,xc]=dmc(N+1,2,L/2);
options = odeset('abstol', 1e-6,'reltol', 1e-4,'maxstep',.6,'stats', 'on');
[t,y]=ode15s(@sbf_sp,[0,2500],y,options,D2,xc);
time=toc

% Show results from figure 8 in the paper
% Plot for a location inside the stimulus region
% for all variables except V_d which is shown 2/3 of
% the way down the cable

% Plot V_sh
figure;
plot(t,y(:,2+2));

% Plot V_d
figure;
plot(t,y(:,N+22));

% Plot Ca
figure;
plot(t,y(:,2*N+2));

% Plot Rss (Scaled)
figure;
plot(t,y(:,3*N+2)/10^6);

% Plot nbar 
figure;
plot(t,y(:,4*N+2));

save figure8  y t;



