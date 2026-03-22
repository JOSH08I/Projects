%% Definition of parameters

% Symbolic variable 
syms s;

% Define the physical parameters
EYoung1 = 120*10^9;
EYoung2 = 120*10^9;
b = 10*10^-3;
h1s = 0.2*10^-3; %thickness of beam 1
h2s = 0.2*10^-3; %thickness of beam 2
l1 = 0.06;       %length of beam 1
l2 = 0.06;       %length of beam 2
rho1s = 8960;    %density of beam 1
rho2s = 8960;  % density of beam 2
M = 0.3;               %primary system mass
m1 = 0.043;               %beam 1 endmass
m2 = 0.040;             %beam 2 endmass
rad = 0.01;              %radius of both endmasses
lm1=m1 / (rho1s*pi*rad^2);%length of beam 1 endmass
offset1=10.5*10^-3;
J1=m1*(lm1^2/3 + rad^2/4 + offset1^2);%MI of beam 1 endmass
lm2=m2 / (rho2s*pi*rad^2); % length of beam 2 endmass
J2=m2*(lm2^2/3 + rad^2/4 + offset1^2); % MI of beam 2 endmass

%Wake Oscillator model parameters
rhof=1.225;
Qbar=2;
Cl0=0.3842;
Cd=1.1856;
St=0.2;
Cm=1;
L=0.3;
d=0.1;
mf=Cm*rhof*pi*d^2/4;

%Piezoelectric parameters PZT - 5H
d31= - 274*10^-12; 
eps33s=25.5*10^-9;
hp = 0.2*10^-3;
bp = b;
Ep = 17*10^9;
rhop = 7500;

% %Piezoelectric parameters - PVDF
% d31= - 30*10^- 12;
% eps33s=13*8.85*10^- 12;
% hp = 0.2*10^- 3;
% bp = 5*10^- 3;
% Ep = 2.5*10^9;
% rhop = 1780*bp*hp;


%calculated parameters
EI1 = EYoung1*b*h1s^3/12  +  (2/3)*bp*Ep*((hp + h1s/2)^3 - h1s^3/8); %%EI of beam 1
EI2 = EYoung2*b*h2s^3/12  +  (2/3)*bp*Ep*((hp + h2s/2)^3 - h2s^3/8); %%EI of beam 2
rho1 = rho1s*b*h1s + 2*rhop*bp*hp ; %linear mass density of beam 1
rho2 = rho2s*b*h2s + 2*rhop*bp*hp; %linear mass density of beam 2

%Calculating Natural Frequencies

syms t1 t2 beta1 beta2;

det_1 = J1 * m1 * beta1^2 * ((EI1 * beta1^4) / rho1)^2 * (2  -  2 * cos(l1 * beta1) * cosh(l1 * beta1))  +  ...
    (EI1^2 * beta1^6 * (2 * rho1  +  cosh(l1 * beta1) * (2 * rho1 * cos(l1 * beta1)  +  ...
    beta1 * ( - 2 * m1  -  2 * J1 * beta1^2) * sin(l1 * beta1))  +  beta1 * (2 * m1  -  2 * J1 * beta1^2) * cos(l1 * beta1) * sinh(l1 * beta1))) / rho1;

det_2 = J2 * m2 * beta2^2 * ((EI2 * beta2^4) / rho2)^2 * (2  -  2 * cos(l2 * beta2) * cosh(l2 * beta2))  +  ...
    (EI2^2 * beta2^6 * (2 * rho2  +  cosh(l2 * beta2) * (2 * rho2 * cos(l2 * beta2)  +  ...
    beta2 * ( - 2 * m2  -  2 * J2 * beta2^2) * sin(l2 * beta2))  +  beta2 * (2 * m2  -  2 * J2 * beta2^2) * cos(l2 * beta2) * sinh(l2 * beta2))) / rho2;

BetaSolutions1 = [];
BetaSolutions2 = [];
rangeStart = 0.1;  % Start at 0.1 to avoid zero solution
rangeEnd = 20;    % Specify the upper limit for Beta search
step = 2;          % Step size for finding different solutions

for i = rangeStart:step:rangeEnd
    sol1 = vpasolve(det_1, beta1, [i, i  +  step]);
    if ~isempty(sol1)
        BetaSolutions1 = [BetaSolutions1; sol1];  % Collect solutions
    end
end
for i = rangeStart:step:rangeEnd
    sol2 = vpasolve(det_2, beta2, [i, i  +  step]);
    if ~isempty(sol2)
        BetaSolutions2 = [BetaSolutions2; sol2];  % Collect solutions
    end
end
omega_beam1 = (double(BetaSolutions1(1))^4 * EI1 / rho1)^0.5; %Natural frequency of the first beam
omega_beam2 = (double(BetaSolutions2(1))^4 * EI2 / rho2)^0.5; %Natural frequency of the second beam

%% Defining Additional Parameters
masstot=M + m1 + m2; %total mass
omega1=50;
k=omega1^2*masstot;
zeta1=0.05;
zeta2=0.05; 
c1 = zeta1*2*masstot*omega1;
c21 =  zeta2*2*m1*omega_beam1;
c22 =  zeta2*2*m2*omega_beam2;
alpha = omega1 / omega_beam1;
beta = omega1 / omega_beam2;

% Mode Shape calculation
ortho1 = (EI1 * t1^2 * beta1^3 * (J1 * beta1^3 * (J1 * l1 * beta1^4  -  2 * rho1)  +  ...
        l1 * beta1 * ( - J1^2 * beta1^6  +  (0  +  1i) * rho1^2) * cos((1  +  1i) * l1 * beta1)  +  ...
        beta1 * (0.5 * J1^2 * l1 * beta1^6  +  J1 * beta1^2 * rho1  -  0.5 * l1 * rho1^2) * cos(2 * l1 * beta1)  +  ...
        l1 * beta1 * ( - J1^2 * beta1^6  -  (0  +  1i) * rho1^2) * cosh((1  +  1i) * l1 * beta1)  +  ...
        (0.5 * J1^2 * l1 * beta1^7  +  J1 * beta1^3 * rho1  +  0.5 * l1 * beta1 * rho1^2) * cosh(2 * l1 * beta1)  +  ...
        (0.5  -  0.5i) * J1^2 * beta1^6 * sin((1  +  1i) * l1 * beta1)  -  (1  +  1i) * J1 * l1 * beta1^4 * rho1 * sin((1  +  1i) * l1 * beta1)  +  ...
        (0.5  +  0.5i) * rho1^2 * sin((1  +  1i) * l1 * beta1)  -  0.25 * J1^2 * beta1^6 * sin(2 * l1 * beta1)  +  ...
        J1 * l1 * beta1^4 * rho1 * sin(2 * l1 * beta1)  +  0.25 * rho1^2 * sin(2 * l1 * beta1)  -  ...
        (0.125  -  0.125i) * J1^2 * beta1^6 * sin((2  +  2i) * l1 * beta1)  +  ...
        (0.125  +  0.125i) * rho1^2 * sin((2  +  2i) * l1 * beta1)  +  (0.5  -  0.5i) * J1^2 * beta1^6 * sinh((1  +  1i) * l1 * beta1)  +  ...
        (1  +  1i) * J1 * l1 * beta1^4 * rho1 * sinh((1  +  1i) * l1 * beta1)  -  (0.5  +  0.5i) * rho1^2 * sinh((1  +  1i) * l1 * beta1)  -  ...
        0.25 * J1^2 * beta1^6 * sinh(2 * l1 * beta1)  -  J1 * l1 * beta1^4 * rho1 * sinh(2 * l1 * beta1)  -  ...
        0.25 * rho1^2 * sinh(2 * l1 * beta1)  -  (0.125  -  0.125i) * J1^2 * beta1^6 * sinh((2  +  2i) * l1 * beta1)  -  ...
        (0.125  +  0.125i) * rho1^2 * sinh((2  +  2i) * l1 * beta1))) / ...
        (rho1 * cos(l1 * beta1)  +  rho1 * cosh(l1 * beta1)  -  J1 * beta1^3 * sin(l1 * beta1)  -  J1 * beta1^3 * sinh(l1 * beta1))^2 == omega_beam1^2;

ortho2 = (EI2 * t2^2 * beta2^3 * (J2 * beta2^3 * (J2 * l2 * beta2^4  -  2 * rho2)  +  ...
        l2 * beta2 * ( - J2^2 * beta2^6  +  (0  +  1i) * rho2^2) * cos((1  +  1i) * l2 * beta2)  +  ...
        beta2 * (0.5 * J2^2 * l2 * beta2^6  +  J2 * beta2^2 * rho2  -  0.5 * l2 * rho2^2) * cos(2 * l2 * beta2)  +  ...
        l2 * beta2 * ( - J2^2 * beta2^6  -  (0  +  1i) * rho2^2) * cosh((1  +  1i) * l2 * beta2)  +  ...
        (0.5 * J2^2 * l2 * beta2^7  +  J2 * beta2^3 * rho2  +  0.5 * l2 * beta2 * rho2^2) * cosh(2 * l2 * beta2)  +  ...
        (0.5  -  0.5i) * J2^2 * beta2^6 * sin((1  +  1i) * l2 * beta2)  -  (1  +  1i) * J2 * l2 * beta2^4 * rho2 * sin((1  +  1i) * l2 * beta2)  +  ...
        (0.5  +  0.5i) * rho2^2 * sin((1  +  1i) * l2 * beta2)  -  0.25 * J2^2 * beta2^6 * sin(2 * l2 * beta2)  +  ...
        J2 * l2 * beta2^4 * rho2 * sin(2 * l2 * beta2)  +  0.25 * rho2^2 * sin(2 * l2 * beta2)  -  ...
        (0.125  -  0.125i) * J2^2 * beta2^6 * sin((2  +  2i) * l2 * beta2)  +  ...
        (0.125  +  0.125i) * rho2^2 * sin((2  +  2i) * l2 * beta2)  +  (0.5  -  0.5i) * J2^2 * beta2^6 * sinh((1  +  1i) * l2 * beta2)  +  ...
        (1  +  1i) * J2 * l2 * beta2^4 * rho2 * sinh((1  +  1i) * l2 * beta2)  -  (0.5  +  0.5i) * rho2^2 * sinh((1  +  1i) * l2 * beta2)  -  ...
        0.25 * J2^2 * beta2^6 * sinh(2 * l2 * beta2)  -  J2 * l2 * beta2^4 * rho2 * sinh(2 * l2 * beta2)  -  ...
        0.25 * rho2^2 * sinh(2 * l2 * beta2)  -  (0.125  -  0.125i) * J2^2 * beta2^6 * sinh((2  +  2i) * l2 * beta2)  -  ...
        (0.125  +  0.125i) * rho2^2 * sinh((2  +  2i) * l2 * beta2))) / ...
        (rho2 * cos(l2 * beta2)  +  rho2 * cosh(l2 * beta2)  -  J2 * beta2^3 * sin(l2 * beta2)  -  J2 * beta2^3 * sinh(l2 * beta2))^2 == omega_beam2^2;

orthonum1 = subs(ortho1, beta1, BetaSolutions1(1));
orthonum2 = subs(ortho2, beta2, BetaSolutions2(1));

t1num=abs(vpasolve(orthonum1,t1));
t2num=abs(vpasolve(orthonum2,t2));

A1num=t1num(1);

B1=t1num(1)*( - 1 * J1 * beta1^3 * cos(l1 * beta1)  +  1 * J1 * beta1^3 * cosh(l1 * beta1) ...
     -  1 * rho1 * sin(l1 * beta1)  -  1 * rho1 * sinh(l1 * beta1)) / ...
    (1 * rho1 * cos(l1 * beta1)  +  1 * rho1 * cosh(l1 * beta1) ...
     -  1 * J1 * beta1^3 * sin(l1 * beta1)  -  1 * J1 * beta1^3 * sinh(l1 * beta1));
B1num=subs(B1, beta1, BetaSolutions1(1));

C1num= - A1num;
D1num= - B1num;

A2num=t2num(1);
B2 = t2num(1)*( - 1 * J2 * beta2^3 * cos(l2 * beta2)  +  1 * J2 * beta2^3 * cosh(l2 * beta2) ...
     -  1 * rho2 * sin(l2 * beta2)  -  1 * rho2 * sinh(l2 * beta2)) / ...
    (1 * rho2 * cos(l2 * beta2)  +  1 * rho2 * cosh(l2 * beta2) ...
     -  1 * J2 * beta2^3 * sin(l2 * beta2)  -  1 * J2 * beta2^3 * sinh(l2 * beta2));
B2num=subs(B2, beta2, BetaSolutions2(1));
C2num= - A2num;
D2num= - B2num;

beta1num=BetaSolutions1(1);
beta2num=BetaSolutions2(1);

phi1 = A1num*sin(beta1num*s) + B1num*cos(beta1num*s) + C1num*sinh(beta1num*s) + D1num*cosh(beta1num*s);
phi2 = A2num*sin(beta2num*s) + B2num*cos(beta2num*s) + C2num*sinh(beta2num*s) + D2num*cosh(beta2num*s);
%% Calculation of  parameters
% Precompute integrals and numerical values
eta11 = double(int(c21*phi1^2, s, 0, l1));
eta12 = double(EI1*(int(phi1*diff(phi1,s)^2*diff(phi1,s,4) + 4*phi1*diff(phi1,s)*diff(phi1,s,2)*diff(phi1,s,3) + phi1*diff(phi1,s,2)^3,s,0,l1)));
eta13 = double(m1 * (int(diff(phi1, s)^2, s, 0, l1))^2);
eta14 = double(m1 * int(diff(phi1, s)^2, s, 0, l1));
eta21 = double(int(c22*phi2^2, s, 0, l2));
eta22 = double(EI2*(int(phi2*diff(phi2,s)^2*diff(phi2,s,4) + 4*phi2*diff(phi2,s)*diff(phi2,s,2)*diff(phi2,s,3) + phi2*diff(phi2,s,2)^3,s,0,l2)));
eta23 = double(m2 * (int(diff(phi2, s)^2, s, 0, l2))^2);
eta24 = double(m2 * int(diff(phi2, s)^2, s, 0, l2));

X0 = 0.001;
epsilon1 = eta14 * X0;
epsilon2 = eta24 * X0;
R1 = (masstot) / m1;
R2 = (masstot) / m2;
mu21 = eta11 / (2 * epsilon1);
mu22 = eta21 / (2 * epsilon2);
mu1 = c1 / (2 * epsilon1 * (masstot));
delta1 = (masstot) * eta12 / eta14^2;
delta2 = (masstot) * eta22 / eta24^2;
%% Forcing and piezo params
vp1= - bp*Ep*d31*(hp + h1s);
vp2= - bp*Ep*d31*(hp + h2s);
thetap1=double(vp1*(subs(diff(phi1,s),l1)));
thetap2=double(vp2*(subs(diff(phi2,s),l2)));
Cp1=2*eps33s*bp*l1/hp;
Cp2=2*eps33s*bp*l2/hp;
R=10^4;
% %% single result calculation
% Ur=8;
% Uinf=Ur*omega1*d/(2*pi);
% omegaf=2*pi*Uinf*St/d;
% 
% if Ur<6.5
%     epsx=0.05;
%     Ax=4;
% else
%     epsx=0.7;
%     Ax=12;
% end
% tic
% 
%     initial_conditions = [0, 10^- 6, 10^- 6, 0, 0, 0, 0, 0, 0.01, 0];
% 
%     % Time span
%     tspan = [0 500];
% 
%     % Solve the system of ODEs 
%     [t, Y] = ode45(@(t, Y) coupledODEs(t, Y, epsilon1, epsilon2, R1, R2, mu1, mu21, mu22, delta1, delta2, omega_beam1, omega_beam2, omega1, thetap1, thetap2, Cp1, Cp2, R,masstot,rhof,Uinf,L,d,Cl0,Qbar,Cd,X0, Ax, epsx, omegaf), tspan, initial_conditions);
% 
%     % Extract solutions
%     x = Y(:,1);
%     q1 = Y(:,2);
%     q2 = Y(:,3);
%     V1 = Y(:,7);
%     V2 = Y(:,8);
%     Q = Y(:,9);
% % Plot each variable in a separate subplot
% toc;
% %%
% figure;
% plot(t, x);
% xlabel('Time');
% ylabel('x');
% title('Displacement x');
% 
% figure;
% plot(t, q1);
% xlabel('Time');
% ylabel('q1');
% title('Displacement q1');
% 
% figure;
% plot(t, q2);
% xlabel('Time');
% ylabel('q2');
% title('Displacement q2');
% 
% figure;
% plot(t, V1);
% xlabel('Time');
% ylabel('V1');
% title('Voltage V1');
% 
% figure;
% plot(t, V2);
% xlabel('Time');
% ylabel('V2');
% title('Voltage V2');
% 
% figure;
% plot(t, Q);
% xlabel('Time');
% ylabel('Q');
% title('Wake Oscillator Parameter Q');
%% Velocity Sweep 
%Initialize arrays for results
U_values = 1:0.01:10; % Define the range of r values
num_U_values = length(U_values);
max_x_values = zeros(num_U_values, 1);
max_x_dot_values = zeros(num_U_values, 1);
max_q1_values = zeros(num_U_values, 1);
max_q2_values = zeros(num_U_values, 1);
max_V1_values = zeros(num_U_values, 1);
max_V2_values = zeros(num_U_values, 1);
max_Q_values = zeros(num_U_values, 1);
tic
% Loop over the range of r values
for i = 1:num_U_values
    Ur = U_values(i);
    Uinf = Ur*(omega1*d)/(2*pi);
    omegaf=2*pi*Uinf*St/d;

if Ur<6.5
    epsx=0.05;
    Ax=4;
else
    epsx=0.7;
    Ax=12;
end

    % Define time span and initial conditions
initial_conditions = [10^- 6, 10^- 6, 10^-6, 0, 0, 0, 0, 0, 0.001, 0];

    % Time span
    tspan = [0 500];

    % Solve the system of ODEs 
    [t, Y] = ode15s(@(t, Y) coupledODEs(t, Y, epsilon1, epsilon2, R1, R2, mu1, mu21, mu22, delta1, delta2, omega_beam1, omega_beam2, omega1, thetap1, thetap2, Cp1, Cp2, R,masstot,rhof,Uinf,L,d,Cl0,Qbar,Cd,X0, Ax, epsx, omegaf), tspan, initial_conditions);

    % Extract solutions
    x = Y(:,1);
    q1 = Y(:,2);
    q2 = Y(:,3);
    V1 = Y(:,7);
    V2 = Y(:,8);
    Q = Y(:,9);

    % Step 1: Find indices for the time range [2300, 2400]
    time_indices = find(t >= 450 & t <= 500);

    % Step 2: Extract the corresponding values of x, q1, and q2
    x_values = x(time_indices);
    x_dot_values = x(time_indices);
    q1_values = q1(time_indices);
    q2_values = q2(time_indices);
    V1_values = V1(time_indices);
    V2_values = V2(time_indices);
    Q_values = Q(time_indices);

    % Step 3: Compute the maximum values
    max_x_values(i) = max(x_values);
    max_x_dot_values(i) = max(x_dot_values);
    max_q1_values(i) = max(q1_values);
    max_q2_values(i) = max(q2_values);
    max_V1_values(i) = max(V1_values);
    max_V2_values(i) = max(V2_values);
    max_Q_values(i) = max(Q_values);
end
toc
%%
x_dim=max_x_values*X0;
x_dot_dim=max_x_dot_values*X0;
Q_dim=max_Q_values*X0;
y1_dim=max_q1_values*((M + m1 + m2)^0.5)*X0*subs(phi1,s,l1);
y2_dim=max_q2_values*((M + m1 + m2)^0.5)*X0*subs(phi2,s,l2);
V1_dim=max_V1_values*((M + m1 + m2)^0.5)*X0;
V2_dim=max_V2_values*((M + m1 + m2)^0.5)*X0;
Uinfvalues=U_values*(omega1*d)/(2*pi);
Uinfvalues=Uinfvalues';

fig = figure;
sSize = get(0, 'ScreenSize');
set(fig, 'Position', [1, 1, min(sSize(3:4)), min(sSize(3:4))]);
hold on;
plot(Uinfvalues, x_dim , 'LineWidth', 1.5);
plot(Uinfvalues, y1_dim, 'LineWidth', 1.5);
plot(Uinfvalues, y2_dim, 'LineWidth', 1.5);
hold off;
% Set axis labels and legend with LaTeX interpreter
xlabel('$$U_{\infty}$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$Amplitude~(m)$', 'FontSize', 20, 'Interpreter', 'latex');

% Add a legend with LaTeX font
legend({'$Primary$','$Beam~1$','$Beam~2$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');

% Set tick label font to LaTeX
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';  % Ensure LaTeX font is used for tick labels
ylim([0 0.008]);

fig = figure;
sSize = get(0, 'ScreenSize');
set(fig, 'Position', [1, 1, min(sSize(3:4)), min(sSize(3:4))]);
hold on;
plot(Uinfvalues, V1_dim, 'LineWidth', 1.5);
plot(Uinfvalues, V2_dim, 'LineWidth', 1.5);

% Set axis labels and legend with LaTeX interpreter
xlabel('$U_{\infty}$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$Voltage~(V)$', 'FontSize', 20, 'Interpreter', 'latex');

% Add a legend with LaTeX font
legend({'$Beam~1$','$Beam~2$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');

% Set tick label font to LaTeX
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';  % Ensure LaTeX font is used for tick labels

Power1=V1_dim.^2./R;
Power2=V2_dim.^2./R;

fig = figure;
sSize = get(0, 'ScreenSize');
set(fig, 'Position', [1, 1, min(sSize(3:4)), min(sSize(3:4))]);
hold on;
plot(Uinfvalues, Power1, 'LineWidth', 1.5);
plot(Uinfvalues, Power2, 'LineWidth', 1.5);

% Set axis labels and legend with LaTeX interpreter
xlabel('$$U_{\infty}$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$Power~(W)$', 'FontSize', 20, 'Interpreter', 'latex');

% Add a legend with LaTeX font
legend({'$Beam~1$','$Beam~2$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');
% Set tick label font to LaTeX
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';  % Ensure LaTeX font is used for tick labels
%%
Uinfvalues=U_values*(omega1*d)/(2*pi);
Uinfvalues=Uinfvalues';
omegaf=2*pi*Uinfvalues*St/d;

F = (1 / (2 * masstot)) * rhof .* Uinfvalues.^2 * L * d .* ...
    (Q_dim * Cl0 / Qbar - (Cd .* x_dot_dim ./ Uinfvalues)) .* ...
    sqrt(1 + (x_dot_dim ./ Uinfvalues).^2);
fig = figure;
sSize = get(0, 'ScreenSize');
set(fig, 'Position', [1, 1, min(sSize(3:4)), min(sSize(3:4))]);
hold on;
plot(Uinfvalues, F,"LineWidth",1.5);
% Set axis labels and legend with LaTeX interpreter
xlabel('$U_r~(m/s)$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$Force~(N)$', 'FontSize', 20, 'Interpreter', 'latex');

% Add a legend with LaTeX font
% legend({'$Beam~1$','$Beam~2$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best');

% Set tick label font to LaTeX
ax = gca;
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';  % Ensure LaTeX font is used for tick labels

%plot(U_values,Uinfvalues,"Linewidth",1.5)
%%
% Define the system of ODEs for numerical solution
function dqdt = coupledODEs(t, Y, epsilon1, epsilon2, R1, R2, mu1, mu21, mu22, delta1, delta2, omegabeam1, omegabeam2, omega1, thetap1, thetap2, Cp1, Cp2, R, masstot,rhof,Uinf,L,d,Cl0,Qbar,Cd,X0, Ax, epsx , omegaf)
    % Unpack variables
    x = Y(1);
    q1 = Y(2);
    q2 = Y(3);
    x_dot = Y(4);
    q1_dot = Y(5);
    q2_dot = Y(6);
    V1 = Y(7);
    V2 = Y(8);
    Q=Y(9);
    Q_dot=Y(10);

    % Define ODE system
    x_ddot = (1/(2*(masstot))*rhof*Uinf^2*L*d*(Q*Cl0/Qbar - Cd*x_dot/Uinf)*(1 + (x_dot*X0/Uinf)^2)^(1/2) - 2 * epsilon1 * mu1 * x_dot  -  omega1^2 * x  +  epsilon1 * (q1_dot^2) +  epsilon2 * (q2_dot^2) ...
               - epsilon1 * q1 * ( 2 * epsilon1 * mu21 * q1_dot  +  omegabeam1^2 * q1  +  delta1 * epsilon1^2 * q1^3  +  R1 * epsilon1^2 * q1 * q1_dot^2 + thetap1*V1)/(1  +  R1 * epsilon1^2 * q1^2)...
               - epsilon2 * q2 * ( 2 * epsilon2 * mu22 * q2_dot  +  omegabeam2^2 * q2  +  delta2 * epsilon2^2 * q2^3  +  R2 * epsilon2^2 * q2 * q2_dot^2 + thetap2*V2)/(1  +  R2 * epsilon2^2 * q2^2))...
               /(1 - (epsilon1^2*q1^2)/(1  +  R1 * epsilon1^2 * q1^2) - (epsilon2^2*q2^2)/(1  +  R2 * epsilon2^2 * q2^2));
     
    q1_ddot = ( - 2 * epsilon1 * mu21 * q1_dot  -  omegabeam1^2 * q1  -  delta1 * epsilon1^2 * q1^3  -  R1 * epsilon1^2 * q1 * q1_dot^2  -  (thetap1*V1)  +  epsilon1 * q1 * x_ddot) / (1  +  R1 * epsilon1^2 * q1^2);
    
    q2_ddot = ( - 2 * epsilon2 * mu22 * q2_dot  -  omegabeam2^2 * q2  -  delta2 * epsilon2^2 * q2^3  -  R2 * epsilon2^2 * q2 * q2_dot^2  -  (thetap2*V2)   + epsilon2 * q2 * x_ddot)/ (1  +  R2 * epsilon2^2 * q2^2);

    V1dot = thetap1*q1_dot/Cp1 - V1/(R*Cp1);

    V2dot = thetap2*q2_dot/Cp2 - V2/(R*Cp2);

    Q_ddot = (Ax/d)*x_ddot - omegaf*epsx*(Q^2*X0^2-1)*Q_dot - omegaf^2*Q;


    % Return derivatives
    dqdt = [x_dot; q1_dot; q2_dot; x_ddot; q1_ddot; q2_ddot; V1dot; V2dot; Q_dot; Q_ddot];
end




