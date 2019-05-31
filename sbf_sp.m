% In this program we re-write the system of PDE's in terms of ODE's by 
% using the spectral method. We estimate the second derivative term 
% in the the first equation for V_d using the spectral method. 


function yp=sbf_sp(t,y,D2,xc)
sy=length(y);                   % Length of the initial vector
N=sy/8;                         % number of compartments
L=3;                            

        % Define vectors Vd, Vsh, Ca, Rss, nbar, n, m, and h
Vd=  y(1:N);                    % Dendrite Voltage
Vsh= y(N+1:2*N);                % Spine head Voltage
Ca=  y(2*N+1:3*N);              % Calcium concentration
Rss= y(3*N+1:4*N);              % Resistance of spine stem
nbar=y(4*N+1:5*N);              % Spine Density

        % n, m and n are the variable used in the Hodgkin huxley model.
n=   y(5*N+1:6*N);              
m=   y(6*N+1:7*N);
h=   y(7*N+1:8*N); 

                    % Parameter List
                    
A_sh = 1.31e-8;      % Surface area of each spine head cm^2
C_crit = 300;        % Critical intraspine Calcium level (nM)
C_m = 10^(-3);       % Specific membrane Capacitance (mF/cm^2)
C_min = 5;           % Calcium Lower bound (nM)
d = 0.000036;        % Dendritic Cable diameter (cm)
eps_1 = 3e-3;        % Rate of Change in Ca equation (ms^-1)
eps_2 = 7.5e-5;      % Rate of Change in Rss equation (ms^-1)
eps_3 = 1e-5;        % Rate change in nbar equation
gamma = 2.5;         % Channel Density scale factor 
L = 3;               % Dimensionless length of the cable
gbar_Na = .120;      % Maximal Sodium conductance (S/cm^2)
gbar_K = .036;       % Maximal potassium concuctance (S/cm^2);
gbar_L = .0003;      % Maximal Leackage Conductance (S/cm^2)

kappa_c = 1e-9;      % Scale Factor (Calcium Model), (mA*ms/nM)

R_i = 70;            % Specific cytoplasmic resistivity (Ohm-cm)
R_m = 2500;          % Passive membrane resistence (Ohm-cm^2)
R_max = 1000000000;  % Stem resistence upper bound (Ohm)
R_min = 30000000;    % Stem resistence lower bound (Ohm)
R_sh = 1.02e11;      % Resistence of each spine head (Ohm)
V_Na = 115;          % Sodium reversal potential (mV)
V_K = -12;           % Potassium reversal potential (mV)
V_L = 10.5989;       % Leakage reversal potential (mV)
V_syn = 100;         % synaptic reversal potential (mV)

temp=22;             % Temperature parameters
phi=3^((temp-6.3)/10);

nbar_max = 100;     % Spine density upper bound (5-6 spine/10 micro m)
nbar_min = 16;      % Spine density lower bound
Ca_1 = 30;          % Lower bound of Ca concentration where LTD changes to LTP
Ca_2 = 300;         % Upper bound of Ca concentration 

% Define I1 and I2 appearing in boundary conditions
I_1 = 0;            % Injected current in the dendrite
I_2 = 0;            % Released current on the opposite side of dendrite

% Define tau_m, lambda, R_inf, and C_sh
tau_m=R_m*C_m;                          % Membrane time constant (ms)
lambda=sqrt((R_m*d)/(4*R_i));           % Space constant (cm)
R_inf=R_m/(pi*lambda*d);                % Specific input resistance (Ohm)
C_sh=A_sh*C_m;                          % Compartment specific capacitance (mF)

% Initialize vector yp
yp=zeros(8*N,1);

% Define yp(1:N): the equation for V_d  (using spectral method)
x=xc(2:N+1);
dx=xc(2)-xc(1);
D2n=D2(2:N+1,2:N+1);
D2n(:,1)=D2n(:,1)+D2(2:N+1,1);
D2n(:,N)=D2n(:,N)+D2(2:N+1,N+2);
C=dx*R_inf*(I_1*D2(2:N+1,1)+I_2*D2(2:N+1,N+2));
Iss=(Vsh-Vd)./Rss;
yp(1:N)=D2n*Vd+R_inf*nbar.*Iss+C;
yp(1:N)=yp(1:N)/tau_m;

% Function to determine the stimulus along the cable

vsyn=zeros(N,1);
for k=1:N
   vsyn(k)=gsyn(x(k),t);  
end
Isyn=vsyn.*(Vsh-V_syn);

% Define yp(N+1:2*N): the equation for V_sh.

Iion=gamma*A_sh*(gbar_Na*(Vsh-V_Na).*(m.^3).*h + ...
            gbar_K*(Vsh-V_K).*(n.^4)+gbar_L*(Vsh-V_L));
yp(N+1:2*N)=(-Iion-Isyn-Iss)/C_sh;

% Define yp(2*N+1:3*N): the equation for Ca.

yp(2*N+1:3*N)=-eps_1*(Ca-C_min)+abs(Iss)/kappa_c;

% Define yp(3*N+1:4*N): the equation for Rss.

yp(3*N+1:4*N)=-eps_2*(Rss-R_min).*(1-Rss/R_max).* ...
                     (Ca/Ca_1 - 1).*(Ca/Ca_2-1).*(Ca/C_min - 1);

% Define yp(4*N+1:5*N): the equation for nbar.

yp(4*N+1:5*N) =-eps_3*(Ca/Ca_1-1).*(Ca/Ca_2-1).*(Ca/C_min - 1).*(1-nbar/nbar_max).*(nbar-nbar_min);

% Define yp(5*N+1:6*N): the equation for n.

alfn=phi*0.01*(-Vsh+10)./(exp((-Vsh+10)/10)-1);
betn=phi*0.125*exp(-Vsh/80);
yp(5*N+1:6*N)=alfn.*(1-n)-betn.*n;

% Define yp(6*N+1:7*N): the equation for m.

alfm=phi*0.1*(-Vsh+25)./(exp((-Vsh+25)/10)-1);
betm=phi*4*exp(-Vsh/18);
yp(6*N+1:7*N)=alfm.*(1-m)-betm.*m;

% Define yp(7*N+1:8*N): the equation for h.

alfh=phi*0.07*exp(-Vsh/20);
beth=phi*1./(exp((-Vsh+30)/10)+1);
yp(7*N+1:8*N)=alfh.*(1-h)-beth.*h;

