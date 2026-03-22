clear;
clc;
E = input('Enter the value of Modulus of Elasticity of material in Pa: ');
I = input('Enter value of area moment of inertia of beam column: ');
L= input('Enter length of beam column: ');
invalidCondition = true; %flag for the invalid condition
while invalidCondition
    flag_1=true;
    fprintf('First Boundary Condition \n Press 1 for Fixed \n Press 2 for hinged \n Press 3 for guided \n');
    BC1 = input('Enter the first Boundary condition: ');
    fprintf('Second boundary condition: \n Press 1 for Fixed \n Press 2 for hinged \n Press 3 for guided \n Press 4 for Free \n');
    BC2 = input('Enter the second boundary condition: ');

    if (BC1 == 3 || BC1 == 2) && BC2 == 4
        fprintf('Invalid Boundary Condition, Please input again\n');
    elseif BC1==3 && BC2==3
        fprintf('Invalid Boundary Condition, Please input again\n');
    else
        invalidCondition = false; % Set the flag to exit the loop
    end
end
n=input("enter number of terms you want: ");
%**********************Solves the Eigenvalue problem************************
syms A B C D K x;
l=1;
assume(K~=0);
%charachteristic equation
y(x)= A*sin(K*x) + B*cos(K*x) + C*x + D;
%derivatives up to the third one
dy_dx=diff(y,x);
d2y_dx2=diff(y,x,2);
d3y_dx3=diff(y,x,3);
% BCs according to the applied boundary conditions
if BC1==1 && BC2==1
    bc1=y;
    bc2=dy_dx;
    bc3=y;
    bc4=dy_dx;
    ig=2*pi;
elseif  BC1==1 && BC2==2
    bc1=y;
    bc2=dy_dx;
    bc3=y;
    bc4=d2y_dx2;
    ig=4;
elseif  BC1==1 && BC2==3
    bc1=y;
    bc2=dy_dx;
    bc3=dy_dx;
    bc4=d3y_dx3+(K^2)*dy_dx;
    ig=pi;
elseif BC1==1 && BC2==4
    bc1=y;
    bc2=dy_dx;
    bc3=d2y_dx2;
    bc4=d3y_dx3+(K^2)*dy_dx;
    ig=pi/2;
elseif BC1==2 && BC2==1
    bc1=y;
    bc2=d2y_dx2;
    bc3=y;
    bc4=dy_dx;
    ig=2*pi;
elseif BC1==2 && BC2==2
    bc1=y;
    bc2=d2y_dx2;
    bc3=y;
    bc4=d2y_dx2;
    ig=pi;
elseif BC1==2 && BC2==3
    bc1=y;
    bc2=d2y_dx2;
    bc3=dy_dx;
    bc4=d3y_dx3+(K^2)*dy_dx;
    ig=pi/2;
elseif BC1==3 && BC2==1
    bc1=dy_dx;
    bc2=d3y_dx3+(K^2)*dy_dx;
    bc3=y;
    bc4=dy_dx;
    ig=pi;
elseif BC1==3 && BC2==2
    bc1=dy_dx;
    bc2=d3y_dx3+(K^2)*dy_dx;
    bc3=y;
    bc4=d2y_dx2;
    ig=pi/2;
end
%making equations out of the boundary conditions
e1=subs(bc1,0)==0;
e2=subs(bc2,x,0)==0;
e3=subs(bc3,x,l)==0;
e4=subs(bc4,x,l)==0;
%Converting the equations into the form Rx=q
eqns=[e1,e2,e3,e4];
vars=[A B C D];
[R,q]=equationsToMatrix(eqns,vars);
% Display the matrix R
disp('Matrix R:');
disp(R);
%calculate the determinant of R
detR=det(R);
disp('Determinant of R is');
disp(detR);
%convert the matlab equation from symbolic to numeric;
f=matlabFunction(detR);
%find the first n roots of the equation
roots_found = zeros(1, n);
for i = 1:n
    if i==1
        initial_guess=ig;
    else
        initial_guess=roots_found(i-1)+pi;
    end
    % Find the roots one by one starting from an initial guess
    roots_found(i) = fzero(f, initial_guess);
end
disp("The first n positive roots are:");
disp(roots_found);
P_sol=zeros(1,n);
for i=1:n
   P_sol(i)=(roots_found(i))^2*(E*I)/L^2;
end   
%*****************************ENERGY METHOD***************************
sym_vars = sym('a', [1 n]); % Creating an array of symbolic variables a1,a2,a3....an;
syms x P;
phi_x=0;
%declaring admissible functions for respective BCs
if BC1==1 && BC2==1
    for i=1:n
        term=sym_vars(i)*(x^(i+1))*(x-L)^2;
        phi_x=phi_x+term;
    end
elseif  BC1==1 && BC2==2
    for i=1:n
        term=sym_vars(i)*(x^(i+1))*(x-L);
        phi_x=phi_x+term;
    end
elseif  BC1==1 && BC2==3
    for i=1:n
        %term=sym_vars(i)*((x^(i+1))*(x-3*L/2)^i);
        term=sym_vars(i)*(1-cos(pi*x/L))^i;
        phi_x=phi_x+term;
    end
elseif BC1==1 && BC2==4
    for i=1:n
        term=sym_vars(i)*(x^(i+1));
        phi_x=phi_x+term;
    end
elseif BC1==2 && BC2==1
    for i=1:n
        term=sym_vars(i)*(x)*(x-L)^(i+1);
        phi_x=phi_x+term;
    end
elseif BC1==2 && BC2==2
     for i=1:n
        term=sym_vars(i)*(x^(i))*(x-L);
        phi_x=phi_x+term;
    end
elseif BC1==2 && BC2==3
    for i=1:n
        term=sym_vars(i)*((x)*(x-2*L))^i;
        phi_x=phi_x+term;
    end
elseif BC1==3 && BC2==1
    for i=1:n
        term=sym_vars(i)*(cos(pi*x/L)+1)^i;
        phi_x=phi_x+term;
    end
elseif BC1==3 && BC2==2
    for i=1:n
        term=sym_vars(i)*((x-L)*(x+L))^i;
        phi_x=phi_x+term;
    end
end
fprintf("The admissible Function is: \n");
disp(phi_x);

dphi_dx=diff(phi_x,x);
d2phi_dx2=diff(phi_x,x,2);
%Computing potential and work done
U=((E*I)/2)*int((d2phi_dx2)^2,x,0,L);
W=(P/2)*int((dphi_dx)^2,x,0,L);
PI=U-W; %computing total potential
eqns = cell(1, n); 
%Differentiating the solution with respect to each symbolic variable from
%a1 to an
for i = 1:n
    eqns{i} = diff(PI, sym_vars(i)) == 0; 
end
%converting the equations to matrix form
[V,w]=equationsToMatrix(eqns,sym_vars);
det_V=det(V);
%solving the value of Pcr symbolically and converting it to numerical
P_sol_EM=solve(det_V==0,P);
P_num_EM=double(P_sol_EM);
fprintf("The energy method solution \n");
P_num=zeros(1,n);
for i=1:n
    P_num(i)=P_num_EM(i);
end
P_num_sort=sort(P_num);
disp(P_num_sort);
fprintf("The eigenvalue solution \n")
disp(P_sol);
%Calculating error
error=zeros(1,n);
for i=1:n
    error(i)=((P_sol(i)-P_num_sort(i))/(P_num_sort(i)))*100;
end
fprintf("the error percentage \n");
disp(error);
%Plotting the Eigenvalue and energy method solutions to each other
plot((1:n),P_num_sort,'b','LineWidth',1.5);
xlabel('Number of Terms');
ylabel('Pcritical');
hold on;
plot((1:n),P_sol,'r','LineWidth',1.5);
hold off;
legend('Energy Method', 'Eigenvalue Solution')