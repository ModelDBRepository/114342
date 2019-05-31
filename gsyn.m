% This program generates the frequancy stimulation  used in the second
% equation for V_sh in the term I_syn
 
function y=gsyn(x,t)

% Set the period

if(t<=500)
   period = 100;         % 10 Hz
  end
if (t>500 & t<1500) 
   period = 8;           % 125 Hz
end
if (t>=1500)
   period = 100;         % 10 Hz
end

% Other parameters

Dx=0.2;
tp = 0.2;               % time to peak for g_syn (ms)
gp = 0.074*10^(-9);     % Peak synaptic conductance  (S);
tloc=mod(t,period);

% Set the conductance for the region of synaptic input

if (x>=0 & x <= Dx)
   y=gp*(tloc/tp)*exp(1-tloc/tp);
else
   y=0;
end
