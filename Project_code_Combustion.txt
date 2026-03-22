clear;
clc;
rho_g=1.2;
rho_s=2650;
C_g=1100;
C_s=900;
k_g=0.052; %increased to 52 in the last example
k_s=4; %increased to 400 in the last example
h_v=40;
S=2000;
x_s=(0.5)*1e-3;   % 1, 0.5, 0.25 grids checked
x_g=x_s;    
L_s=0.05;
L_g=2;
B=L_s/x_s;   %solid grid points
A=L_g/x_g;  %gas grid points

t=(0.5)*0.0001; 
K=10/t;
Ts=zeros(B,K);
Tg=zeros(A,K);
Tg(:,:)=473;
Ts(:,:)=473;
u=1e-3;%performed for 3 velocities, 1e-3, 1e-6 and 1e-2 for different cases. 
phi=rho_g*C_g;
psi=rho_s*C_s;
for n=1:K
    for m=2:B-1
    Tg(m,n+1)=Tg(m,n)+((k_g*t)/(phi*x_g))*(Tg(m+1,n)-2*Tg(m,n)+Tg(m-1,n))-u*t*(Tg(m,n)-Tg(m-1,n))/x_g-h_v*t*(Tg(m,n)-Ts(m,n))+S*t;
    Ts(m,n+1)=Ts(m,n)+(k_s/psi)*(t/(x_s^2))*(Ts(m+1,n)-2*Ts(m,n)+Ts(m-1,n))+((h_v*t)/psi)*(Tg(m,n)-Ts(m,n));
    end
    for m=B-1:A-1
    Tg(m,n+1)=Tg(m,n)+((k_g*t)/(phi*x_g))*(Tg(m+1,n)-2*Tg(m,n)+Tg(m-1,n))-u*t*(Tg(m,n)-Tg(m-1,n))/x_g;
    end
    %Tg(A,n+1)=Tg(A-1,n+1);
end
x1=0:x_g:L_g-x_g;
x2=0:x_s:L_s-x_s;
plot(x1,Tg(:,K),'LineWidth',2,'Color',[0 0 0]);
figure;
plot(x2,Ts(:,K),'LineWidth',2,'Color',[0 0 0]);

