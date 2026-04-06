%% QG 2 layer code part 2: Evolve initial condition in time
% Originally Created 4-3-2026 by Margaret Gregory
% Based on an Octave Code by Glenn Flierl
%% Define friction and filters that can be used to prevent numerical instabilities

frc=exp(-kappa*dt*wv2-r*dt);

switch(nx) % FILTER-> gets rid of small wavenumbers that build up enstrophy. This is where numerical dispersion occurs, but only on a small subsection of the total range
    case 128
      cphi = 0.69*pi;
    case 256 
      cphi = 0.715*pi;
    case 512
      cphi = 0.735*pi;
    otherwise
      cphi = 0.65*pi;
end

wvx=sqrt((k*dx).^2+(l*dy).^2);  %amplitude of wavenumber in kl space
filtr=exp(-18*(wvx-cphi).^7).*(wvx>cphi)+(wvx<=cphi); %filter depends on amplitude-> radially symetric in kl space
filtr(isnan(filtr))=1;

kmax2=((nx/2-1)*k0x).^2; %Also has a truncation-> cuts out wavenumbers that would allias when going from transform space to real
trunc=(wv2<kmax2);

% Note: both filter and truncation might not be necessary: sometimes they do the same thing

%% Time Evolution

t=0;
tc=0;


ts=[];
stat=[];
stat1=[];
P=[];

% Transform PV in both layers into spectral space
qh1=fft2(q1);
%qh2=fft2(q2); %already have qh2 from QG_2_Layer_Parameters.m

dqh1dt_p=0;
dqh2dt_p=0;
dt0=dt;dt1=0;

cm=parula(100);
Y=y(:,1);
X=x(1,:);

z=0;
while t<=tmax+dt/2
  z=z+1;

  % Define PV in real space
  q1=real(ifft2(qh1.*trunc));
  q2=real(ifft2(qh2.*trunc));

  % Calculate streamfunctions in spectral space
  [ph1,ph2]=invert(qh1,qh2,a11,a12,a21,a22);

  % Transform streamfunctions in real space
  p1=real(ifft2(ph1.*trunc));
  p2=real(ifft2(ph2.*trunc));
  
  % Calculate velocities
  [u1,v1]=caluv(ph1,k,l,trunc);
  [u2,v2]=caluv(ph2,k,l,trunc);
  
  %Time evolution (NONLINEAR Equations)
  dqh1dt=-advect(q1,u1+U1,v1,k,l)-beta1*1i*k.*ph1;
  dqh2dt=-advect(q2,u2+U2,v2,k,l)-fft2(v2.*Qy')+rek*wv2.*ph2;

  %Time evolution (LINEAR Equations)
  % dqh1dt=-beta1*1i*k.*ph1;
  % dqh2dt=-fft2(v2.*Qy')+rek*wv2.*ph2;

  % Update PVs
  qh1=frc.*filtr.*(qh1+dt0*dqh1dt+dt1*dqh1dt_p);
  qh2=frc.*filtr.*(qh2+dt0*dqh2dt+dt1*dqh2dt_p);

  dqh1dt_p=frc.*dqh1dt;
  dqh2dt_p=frc.*dqh2dt;

  if tc==0
    dt0=1.5*dt;dt1=-0.5*dt;
  end
  tc=tc+1;
  t=tc*dt;

  % Plot ff and record statistics every few timesteps
  if(rem(tc,tpl)==0) || t==dt
    ts=[ts,t];

    stat=[stat,min(p2,[],'all')];
    stat1=[stat1,.5*mean(u2.^2+v2.^2,'all')];
    
    P=cat(3,P,p2); % Records lower layer streamfunction 
    % U=cat(3,U,u2);
    % V=cat(3,V,v2);

    ff=[p2]; % Plots p2, can be changed ex. ff=[rs(p1),rs(p2)]; will plot both streamfunctions

    imagesc(ff);axis('xy');title("\psi_{2}:"+sprintf(" t=%g",t),'FontSize',15); %,'equal'
    colorbar
    drawnow;

  end

end
