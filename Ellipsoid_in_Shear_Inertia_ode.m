function [T,X,P,F,St2] = Ellipsoid_in_Shear_Inertia_ode(X0,P0,N,nt,St,B,ka,w,A,d,fluid)
%% Particles in waves
% based on Bakhoday Paskyabi 2015 (2D)
%also uses's jeffery torques + euler's eqn to find particle rotation,
%coupled back to flow using brenner's resistance tensor
%solves equations using matlab ode solvers

%inputs
%X0                     %initial position of particle in [kx,ky,kz]'
%P0                     %initial orientation of particle, unit vector
%N                      %number of periods to simulate
%nt                     %number of points per period to simulate
%St                     %Stokes number for a sphere
%B=1;                   %water to particle density
%ka                     %steepness (Also Fr)
%w                      %aspect ratio
%A                      %wave amplitude
%d                      %diameter of equivalent volume sphere
%fluid                  %logical 1 or 0- if you want to evaluate fluid path

%outputs
%T                      %time vector
%X                      %position and velocity matrix (m, m/s) [x y z; u v w];
%P                      %orientation of particle over time [p1 p2 p3]
%F                      %corresponding fluid particle path [x, z]
%St2                    %effective Stokes number for ellipsoid w equivalent volume of sphere

%make sure they are columns
if isrow(X0); X0=X0'; end;
if isrow(P0); P0=P0'; end;


% Things to do and keep in mind
%Remember, Jeffery's equation only works for neutrally buoyant particles?
%is this actually true? work close enough for close to neutral?
%Also, there are limitations to Maxey-Riley --> breaks down for larger St
%Added mass changes for non-spherical particles.
%no Faxen correction
%No history force
%Coupled particle rotations back to drag and advection


%% constants
g=9.81;                         %acceleration due to gravity, (m/s^2)
nu=1.05*10^-6;                  %kinematicviscosity of seawater - temp specific, now use 20 deg C value
rho=1029;                       %density of seawater (%kg/m^3)
mu=rho*nu;                      %dyanamic viscosity (N s/m^2)

%% Wave parameters using Dispersion  for DEEP WATER WAVES, can use findk if using finite depth
k=ka/A;                         %wave number (1/m)
omega=sqrt(g*k);                %wave frequency (1/s)
T=2*pi/omega;                   %wave period


%% Particle Properties
if ~isempty(St) && isempty(d)
    tau=St/omega;               %particle relaxation time for a sphere (s)
    D=sqrt(tau*18*B*nu);        %Diameter of the particle from a spherical definition of Tau
else
    D=d;                        %use input diameter (m)
    tau=D^2/18/B/nu;            %particle relaxation time for a sphere (s)
end

%by geometry- this has been double checked.
V_sphere=(D/2)^3;               %volume of sphere
b=(V_sphere/w)^(1/3);           %minor ellipsoid radius
a=w*b;                          %major ellipsoid radius
V=4/3*pi*a*b*b;                 %volume of ellipsoid - should be the same as sphere? plus 4/3*pi

[ I,alpha0,beta0,gamma0] = MoI_ellipsoid(a,b,B,rho); %moment of Inertia + torque coefficients

%find effective Stokes number using  ellipsoidal tau, (Zhang et al., 2001) and (Challabotla, Zhao and Andersson 2015)
if w>1
    tau_corr=w*log(w+sqrt(w^2-1))/sqrt(w^2-1); %correction where tau_ell=tau_sphere*tau_corr
elseif w==1
    tau_corr=1;
elseif w<1
    tau_corr=4*(w/2)*(pi-2*atan(w*sqrt(1-w^2)))/(sqrt(1-w^2))*b^2/D^2;
end

St2=omega*tau_corr*tau;         %output Stokes number based on wave frequency

%% Model Setup & Run

%Time
if isempty(nt)
    tspan=[0 N*T];              %let ode solver find best time step
else
    tspan=linspace(0,N*T,N*nt); %specify time step
end

%Initial velocity set to fluid velocity:
[u,~]=stokes(omega,k,A,X0(3));

ws=-2/9*g*(D/2)^2/nu*(1/B-1);        %equilibrium settling velocity for sphere

%Initial position of fluid particle (x,z)), only 2D waves so no y movement
F0=X0([1 3]);

%set INITIAL CONDITIONS to the same as fluid (position of particle (x,z), velocity (u,w))
X0=[X0;u(1);0;ws];

%particle trajectory
[T,X,P] = funWaveParticles3D(tspan,X0,A,k,omega, B, g,D,P0,w,nu,I,alpha0,beta0,gamma0,mu,V,a,b);

%fluid trajectory, 2d, just X and Z
if fluid==1
    options = odeset('RelTol',10^-14);
    [T,F]=ode45(@(t,Y)funWave(t,Y,A,k,omega),T,F0,options);
else
    F=[];
end

end
%% FUNCTIONS
function [TT, X,P] = funWaveParticles3D(ts,X0,A,k,omega, B, g,D,P0,w,nu,I,alpha0,beta0,gamma0,mu,V,a,b)
% get particle trajectory! using simplified maxey riley

[R]=ellipsoid_tensor(w);            %resistance tensor of the ellipsoid
Iinv=I\eye(3);%inv(I);
Y0=[X0; P0; 0; 0; 0];               %initial conditions [x;u;p;dpdt]

 %Added mass coefficient
if w~=1
    %based on Lamb
    Cm=[alpha0/(2-alpha0);beta0/(2-beta0);gamma0/(2-gamma0)];
else
    %A sphere
    Cm=1/2;
end

    options = odeset('RelTol',10^-14);
[TT,Y]=ode15s(@(t,y)funwave_ellipsoids(t,y,B,a,b,k,A,omega,nu,D,V,I,alpha0,gamma0,beta0,mu,R,g,w,Iinv,Cm),ts,Y0,options);

X=Y(:,1:6);
P=Y(:,7:12);
end

function [dy]=funwave_ellipsoids(t,y,B,a,b,k,A,omega,nu,D,V,I,alpha0,gamma0,beta0,mu,R,g,w,Iinv,Cm)
%% translation
X=y(1:3);           %particle position
vp=y(4:6);          %particle velocity
p=y(7:9);           %particle orientation
%p=p/norm(p);
omegan=y(10:12);    %particle orientation time derivative

%[u,DuDt,S,~]=airywaves_deep(k,A,omega,X(1),X(3),t);


[u,dudz]=stokes(omega,k,A,X(3));
DuDt=[0 0 0]';
S=.5*[0 0 dudz;
    0 0 0;
    dudz 0 0];
RR=.5*[0 0 dudz;
    0 0 0;
    -dudz 0 0];
%evaluate particle accelerations due to:

%1. buoyancy force:
term1=(1-B)*[0;0; -g];
%2. pressure gradient:
term2=(B)*DuDt;
%3. drag
u3=u-vp;
Rot=vrrotvec2mat(vrrotvec([1 0 0],p)); %rotation relative to x axis
K=(Rot'*R*Rot);          %rotated ellipsoid resistance tensor
term3=3*pi*nu*B*D/V*K*u3;



% Total particle acceleration: Maxey-Riley
dvpdt=term1+term2+term3;
%added mass effects? for a sphere:
%dvpdt=(dvpdt+B/2*DuDt)/(1+B/2);
%added mass effects? for an ellipsoid (LAMB)
if w~=1
dvpdt=(dvpdt+B*abs(Rot*Cm).*DuDt)./(1+abs(Rot*Cm)*B);
else
    dvpdt=(dvpdt+B*Cm.*DuDt)./(1+Cm*B);
end

dx=[vp;dvpdt];

%% orientation
if w~=1
    pdot=cross(omegan,p);  %explicit derivative of orientation vector
    S=(Rot'*S*Rot);        %rotate fluid strain into frame of particle
    RR=(Rot'*RR*Rot);
    %find Jeffery Torques for an ellipsoid - these are dimensional
    M(1)=32*pi*mu/3/(beta0+gamma0)*(RR(2,3)-omegan(1));
    M(2)=16*pi*mu/3/(b^2*gamma0+a^2*alpha0)*((b^2-a^2)*S(1,3)+(b^2+a^2)*(RR(1,3)-omegan(2)));
    M(3)=16*pi*mu/3/(b^2*beta0+a^2*alpha0)*((a^2-b^2)*S(1,2)+(b^2+a^2)*(RR(1,2)-omegan(3)));
    
    %find new angular velocities from Euler's eqn for rigid body motion
    omegadot=Iinv*(M'-cross(omegan,I*omegan));
    
elseif w==1 %No rotation
    pdot=[0 0 0]';
    omegadot=pdot;
end


dy=[dx; pdot; omegadot];

end

%% constant in time and space, resistance tensor of ellipsoid
function [K] = ellipsoid_tensor(w)
%from Loth 2008

if w<1 %oblate exact
    S=acos(w)/sqrt(1-w^2);
    f_par=4/3*w^(-1/3)*(1-w^2)/(w+(1-2*w^2)*S);
    f_per=8/3*w^(-1/3)*(w^2-1)/(w-(3-2*w^2)*S);
    
elseif w>1% && w<=6 
    %prolate approx
    % f_par= (4/5+w/5)*w^(-1/3);
    % f_per=(3/5+2*w/5)*w^(-1/3);
    
    %prolate exact
    f_par=(4/3)*w^(-1/3)*(1-w^2)/(w-(2*w^2-1)*log(w+sqrt(w^2-1))/sqrt(w^2-1));
    f_per=(8/3)*w^(-1/3)*(w^2-1)/(w+(2*w^2-3)*log(w+sqrt(w^2-1))/sqrt(w^2-1));
    
    % elseif w>6 %needle
    % f_par=2/3*w^(2/3)/(log(2*w)-.5);
    % f_per=2*f_par;
elseif w==1
    f_par=1;
    f_per=1;
end


K=zeros(3,3);
K(1,1)=f_par;  %parallel
K(2,2)=f_per;  %perpendicular
K(3,3)=f_per;  %perpendicular
end

%%
function [ I,alpha0,beta0,gamma0] = MoI_ellipsoid(a,b,B,rho )
%get ellipsoid moment of inertia

I=4/3*rho/B*pi*a*b^2*eye(3,3)/5;
I(1,1)=I(1,1)*2*b^2;
I(2,2)=I(2,2)*(b^2+a^2);
I(3,3)=I(3,3)*(b^2+a^2);


w=a/b;

% if w>1
%     k=sqrt(1-1/w^2);
% %k=(w^2-1)/(w^2+1);
% 
% %integrals done by wolfram alpha!!! thx
% % % alpha0=2*(-k^2+sqrt((k^2-1)/k^4)*k^2*acos(sqrt(1/k^2))+1)/(k^2-1)^2;
% % % gamma0=sqrt(1/k^2)*sqr
% % % beta0=gamma0;
% 
% %k is positive:
% alpha0=2/(k^2-1)-2*asec(k)/(k^2-1)^(3/2);
% gamma0=asec(k)/(k^2-1)^(3/2)+1/(k^2-k^4);
% beta0=gamma0;
% 
% alpha=2*(1/(1-w^2)-w*acos(w)/(1-w^2)^(3/2));
% 
% else
%     k=sqrt(1-w^2);
% end


%from Lamb's HYDRODYNAMICS

%prolate, pg 154
if w>=1
    e=sqrt(1-1/w^2);
    alpha0=2*(1-e^2)/e^3*(1/2*log((1+e)/(1-e))-e);
    beta0=1/e^2-(1-e^2)/2/e^3*log((1+e)/(1-e));
    gamma0=beta0;
elseif w<1 %oblate, pg 701
    e=sqrt(1-w^2);
% %     alpha0=sqrt(1-e^2)/e^3*asin(e)-(1-e^2)/e^2;
% %     beta0=alpha0;
% %     gamma0=2/e^2*(1-sqrt(1-e^2)*asin(e)/e); %<-- TYPO IN THE BOOK
    
 %in the book, gamma0 is along the parallel axis?
 
     gamma0=sqrt(1-e^2)/e^3*asin(e)-(1-e^2)/e^2;
     beta0=gamma0;
     alpha0=2/e^2*(1-sqrt(1-e^2)*asin(e)/e); %<-- TYPO IN THE BOOK

end
    

    
    
end

function [u,dudz]=stokes(omega,k,a,z)

u=omega*k*a^2*exp(2*k*z);
dudz=u*2*k;
u=[u 0 0]';


end


