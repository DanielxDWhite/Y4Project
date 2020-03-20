%----------------------------------------------------------------------%
%----Optical Bloch Equation for a two level atom-----------------------%
%----------------------------------------------------------------------%
clc;
clear all;
close all;
%-----------------------------------------------------------------------%
N=1000; % no of time grid points
t=linspace(0,0.25e-6,N); % time array in micro seconds
dt=t(2)-t(1); % step size in time
%---- initial conditions-------------
rho11(1:N)=1;
rho22(1:N)=0;
rho12(1:N)=0;
rho21(1:N)=0;
Omega_rabi=2*pi*10*1e6; % Rabi frequency 2pi*10 Mhz
delta=0; % Detuning frequency
gamma=0.1*Omega_rabi; % linewidth
function rh11=drho11(t,r12,r21,r22,gam,Omega_rabi,delta)
%   Detailed explanation goes here
rh11=1i*0.5*Omega_rabi'*exp(-1i*delta*t)*r12-1i*0.5*Omega_rabi*exp(1i*delta*t)*r21+2*gam*r22;
end
function rh12=drho12(t,r11,r12,r22,gam,Omega_rabi,delta)
%   Detailed explanation goes here
rh12=1i*0.5*Omega_rabi*exp(1i*delta*t)*(r11-r22)-gam*r12;
end
function rh21=drho21(t,r11,r12,r22,gam,Omega_rabi,delta)
%   Detailed explanation goes here
rh21=-1i*0.5*Omega_rabi'*exp(-1i*delta*t)*(r11-r22)-gam*r12;
end
function rh22=drho22(t,r12,r21,r22,gam,Omega_rabi,delta)
%   Detailed explanation goes here
rh22=-1i*0.5*Omega_rabi'*exp(-1i*delta*t)*r12+1i*0.5*Omega_rabi*exp(1i*delta*t)*r21-2*gam*r22;
end

%----Runge Kutta 4th order ----------------------------------------------%
for i=2:N
    
   k1 = drho22(t(i-1),rho12(i-1),rho21(i-1),rho22(i-1),gamma,Omega_rabi,delta)*dt; 
   l1 = drho11(t(i-1),rho12(i-1),rho21(i-1),rho22(i-1),gamma,Omega_rabi,delta)*dt; 
   m1 = drho12(t(i-1),rho11(i-1),rho12(i-1),rho22(i-1),gamma,Omega_rabi,delta)*dt;
   n1 = drho21(t(i-1),rho11(i-1),rho12(i-1),rho22(i-1),gamma,Omega_rabi,delta)*dt;
   k2 = drho22(t(i-1)+0.5*dt,rho12(i-1)+0.5*m1,rho21(i-1)+0.5*n1,rho22(i-1)+0.5*k1,gamma,Omega_rabi,delta)*dt; 
   l2 = drho11(t(i-1)+0.5*dt,rho12(i-1)+0.5*m1,rho21(i-1)+0.5*n1,rho22(i-1)+0.5*k1,gamma,Omega_rabi,delta)*dt; 
   m2 = drho12(t(i-1)+0.5*dt,rho11(i-1)+0.5*l1,rho12(i-1)+0.5*m1,rho22(i-1)+0.5*k1,gamma,Omega_rabi,delta)*dt; 
   n2 = drho21(t(i-1)+0.5*dt,rho11(i-1)+0.5*l1,rho12(i-1)+0.5*m1,rho22(i-1)+0.5*k1,gamma,Omega_rabi,delta)*dt; 
   
   k3 = drho22(t(i-1)+0.5*dt,rho12(i-1)+0.5*m2,rho21(i-1)+0.5*n2,rho22(i-1)+0.5*k2,gamma,Omega_rabi,delta)*dt; 
   l3 = drho11(t(i-1)+0.5*dt,rho12(i-1)+0.5*m2,rho21(i-1)+0.5*n2,rho22(i-1)+0.5*k2,gamma,Omega_rabi,delta)*dt; 
   m3 = drho12(t(i-1)+0.5*dt,rho11(i-1)+0.5*l2,rho12(i-1)+0.5*m2,rho22(i-1)+0.5*k2,gamma,Omega_rabi,delta)*dt; 
   n3 = drho21(t(i-1)+0.5*dt,rho11(i-1)+0.5*l2,rho12(i-1)+0.5*m2,rho22(i-1)+0.5*k2,gamma,Omega_rabi,delta)*dt; 
   
   k4 = drho22(t(i-1)+dt,rho12(i-1)+m3,rho21(i-1)+n3,rho22(i-1)+k3,gamma,Omega_rabi,delta)*dt; 
   l4 = drho11(t(i-1)+dt,rho12(i-1)+m3,rho21(i-1)+n3,rho22(i-1)+k3,gamma,Omega_rabi,delta)*dt; 
   m4 = drho12(t(i-1)+dt,rho11(i-1)+l3,rho12(i-1)+m3,rho22(i-1)+k3,gamma,Omega_rabi,delta)*dt; 
   n4 = drho21(t(i-1)+dt,rho11(i-1)+l3,rho12(i-1)+m3,rho22(i-1)+k3,gamma,Omega_rabi,delta)*dt; 
  
   
   rho22(i) = rho22(i-1)+((k1+2*(k2+k3)+k4)/6);
   rho11(i) = rho11(i-1)+((l1+2*(l2+l3)+l4)/6);
   rho12(i) = rho12(i-1)+((m1+2*(m2+m3)+m4)/6);
   rho21(i) = rho21(i-1)+((n1+2*(n2+n3)+n4)/6);
end
%--------plotting solutions------------------------------------------------
figure(1)
plot(t,rho22,'linewidth',2)
hold on
plot(t,rho11,'k','linewidth',2)
plot(t,abs(rho12),'r','linewidth',2)
plot(t,abs(rho21),'g','linewidth',2)
xlabel('time','fontSize',14);
ylabel('density matrix elements','fontsize',14)
axis([0 t(N) 0 1.1])
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
axis([0 t(N) 0 1.1])
fh = figure(1);
set(fh, 'color', 'white'); 
