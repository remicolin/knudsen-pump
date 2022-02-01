% program Optimisatiion_Knudsen_V4.2

clear all

%---------------------------------------------
% PARAMETERS TO BE MODIFIED:
CD=1.418; %D1/D2
CL=0.6; %L1/Li
a=250.0e-6; %Pipe radius [m]
Kn0=0.2; %Desired initial Knudsen number []
Tc=300.0; %Temperature cold chamber [K]
Th=400; %Temperature hot chamber [K]

dPi=1e-8; %Itération value for pressure gradient [Pa]
err=1e-16; %Minimum error for mass flow between two stages [m3/s]
n=1000; %number of nodes

%---------------------------------------------
% CONSTANTS AND HYPOTHESEES
%Temperature distributions are assumed linears
m=39.0e-3; %Atomic mass (Argon) [kg]
sigmap=1.018; %Viscous slip coefficient []
sigmaT=1.175; %Thermal slip coefficient []
R=8.31446261815324; %Gas constant [J/(mol.K)]

%---------------------------------------------
% CALCULATION OF k2-CONSTANT  
Tref=273.15; %Reference temperature [K]
muref=211.7e-07; %Reference mean-free-path [m]
omega=0.81; %Free-path coefficient []
alpha=1.40; %HS-model coefficient []
k2=(4*alpha*(7-2*omega)*(5-2*omega))/(5*(alpha+1)*(alpha+2)*(sqrt(2*pi)));

%-----------------------------------------------
% INITIAL PARAMETERS
Li=53.0e-3; %Length of pipe [m]
Tc=300.0; %Temperature cold chamber [K]
Th=400; %Temperature hot chamber [K]
Tin=Tc; %[K]
muin=muref*(Tin/Tref)^omega; %[m]
Pin=k2*muin*sqrt(R/m*Tin)/(Kn0*2*a); %[Pa]
deltaP=0; % Pression between two nodes at point of  bifurcation [Pa]
Tm=Tc+(Th-Tc)/2; %Mean temperature [K]

%-----------------------------------------------
% CALCULATIONS
mdot1=1;
mdot2=0;
while abs(mdot1-mdot2)> err
Tx=Tc+(Th-Tc)*CL;
Tx1=Tx-(Th-Tc)/n;
Tx2=Tx+(Th-Tc)/n;
deltaz=Li/n;

%First stage
a1=a;
rho=Pin*m/Tx1/R;
mu=muref*(Tx1/Tref)^omega;
Kni=k2*mu*sqrt(R*Tx1)/(Pin*2*a1);
termP1=-(pi*a1^4*rho*(1+4*sigmap*Kni)*deltaP)/(8*mu*deltaz);
termT1=(sigmaT*mu*pi*a1^2*(Tx-Tx1))/(Tx1*deltaz);
mdot1=termP1+termT1;

%Second stage
a2=a/CD;
rho=Pin*m/Tx/R;
mu=muref*(Tx/Tref)^omega;
Kni=k2*mu*sqrt(R*Tx)/(Pin*2*a2);
termP2=-(pi*a2^4*rho*(1+4*sigmap*Kni)*deltaP)/(8*mu*deltaz);
termT2=(sigmaT*mu*pi*a2^2*(Tx2-Tx))/(Tx*deltaz);
mdot2=2*(termP2+termT2);

%Itération - convergence of deltaP
if mdot1<mdot2
    deltaP=deltaP-dPi;
else
    deltaP=deltaP+dPi;
end

end

%----------------------------------------
% DISPLAY OF RESULTS
    if mdot1<0
        fprintf('Error: \n   mdot is negative. Try changing the value of "CD" \n \n');
    elseif deltaP<0
        fprintf('Error: \n   deltaP is negative. Try changing the value of "CD" \n \n');
    else
        format shortg
        disp('Mass flow Bifurcation [kg/s]:');
        disp(mdot1);
    end

%----------------------------------------------
%CALCULATION AND DISPLAY OF SINGLE TUBE (CAN BE REMOVED)
%A corresponding, single radius is calculated based on the total volume:
ax=sqrt((a1^2*pi*Li*CL+2*a2^2*pi*Li*(1-CL))/(pi*Li));
mu=muref*(Tm/Tref)^omega;
rho=Pin*m/Tm/R;
termP1x=0;
termT1=(sigmaT*mu*pi*ax^2*(Th-Tc))/(Tm*Li);
mdotx=termP1x+termT1;
disp('Mass flow for single tube [kg/s]:');
disp(mdotx);

% end of program
