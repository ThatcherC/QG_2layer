%% QG 2 layer code part 1: Define domain and initial conditions
% Originally Created 4-3-2026 by Margaret Gregory
% Based on an Octave Code by Glenn Flierl

%% Define domain size and resolution + Time step size

r=0; %friction
rek=0; %friction associated with ekman

L=10; % X-domain width
W=L; % Y-domain width: can be set to not equal L

dt=1/32; % time step
tpl=1/dt;   % dictates steps between plots in actual time
tmax=2000; % time numerics run until


U1=0;   %Background zonal velocity of the flow (Upper Layer)
U2=0;   %Background zonal velocity of the flow (Lower Layer)  

kappa=0;    %viscocity term-> when 0, no explicity viscosity. But there is some numerical dispersion    
nx=4*128;ny=4*128;  %Sets Horizontal resolution
dx=L/nx;    %X-grid spacing
dy=W/ny;    %Y-grid spacing

k0x=2*pi/L; %Now defining wavenumbers based on model resolution
k0y=2*pi/W;
[k,l]=meshgrid([0:nx/2,-nx/2+1:-1]*k0x,[0:ny/2,-ny/2+1:-1]*k0y);    % Wave number grid

[x,y]=meshgrid([1/2:1:nx]/nx*L-L/2,[1/2:1:ny]/ny*W-W/2); 


%% Set beta and internal deformation radii
beta1=0;% On f-plane so no beta in upper layer %Could be set to beta+F1*(U1-U2) with beta=1.728e-3  
%beta2=beta-F2*(U1-U2);  

% Topographic beta is important in bottom layer. Derivative of topography
% shows up in time evolution code
yy=y(:,1);
al=1; %Parameter that controls asymmetry of topography
Qy=sech(yy).^2.*(yy<=0)+sech(al*yy).^2.*(yy>0); %Gives tanh(y) topography.

% del=10; %Ratio of layer depths
% Rd=1; % Rossby Radius of Deformation
F1=0; % F1=0 decouples layers. could be set to 1/Rd^2/(1+del);
F2=1; %could be set to F1*del

%% Define Operators used in time stepping code
wv2=(k.*k+l.*l);    %k^2+l^2
det=wv2.*(wv2+F1+F2);
a11=-(wv2+F2)./det; %A matrix tranforms from q in wavenumber space to psi is wavenumber space
a12=-F1./det;
a21=-F2./det;
a22=-(wv2+F1)./det;

% Set first value to 0 in order to eliminate zero mode
a11(1,1)=0; 
a12(1,1)=0;
a21(1,1)=0;
a22(1,1)=0;

%% Define initial q fields

q1=0*x; % Assume PV in upper layer is zero
%q2=0*x; % Can change this is you start knowing the lower layer PV 

% NOTE: This code is written assuming that you know the 2D streamfunction in
% the lower layer

psi2=.1*sech(x).^2.*sech(y).^2; %You know this. Can specify an analytical form or numerical solution for psi2

psihat2=fft2(psi2); %Tranform psi2 into fourier space
qh2=-wv2.*psihat2-F2*(psihat2); %Calculate the lower layer PV (in spectral space) given psihat2 

% Plot psi2
figure;
p=pcolor(x,y,psi2);
colorbar
set(p, 'EdgeColor','none')


% NOW GO RUN QG_2_Layer_Time_Evolution.mat