clc;
clear;
close all;
%F1
T = 1/10;
K=0.5;
J=0.009;
R=3.4;
r=0.0125;
L=0.47;
d=0.15;
Km=150;
a=0.5*(K^2/(J*R));
Gm=tf([K/(J*R)], [1, K^2/(J*R), 0]);
C0 = tf([Km, Km*a], 1);
Gb = tf([(5/7)*9.81*(r/L)],[1, d, 0]);
Gth=feedback(C0*Gm,1);
%Gth=minreal(Gth,0.0001);
G = Gth*Gb;

%F2

%F3
G1pz =c2d(G, T, 'zoh');

G2pz = minreal(G1pz, 0.0001);

%F4
%Controlable
[num, den] = tfdata(G2pz, 'v');
Ac = [0 1 0; 0 0 1; -den(4) -den(3) -den(2)];
Bc = [0; 0; 1];
Cc = [num(4) num(3) num(2)];
Dc = 0;
sc = ss(Ac, Bc, Cc, Dc, T);
[Vc, Dc] = eig(Ac);
V_1c=inv(Vc);
%Observable
Ao = [0 0 -den(4); 1 0 -den(3); 0 1 -den(2)];
Bo = [num(4);num(3);num(2)];
Co = [0 0 1];
Do = 0;
so = ss(Ao, Bo, Co, Do, T);

[Vo, Do] = eig(Ao);
V_1o=inv(Vo);

%F5


%Controllable form

y_kc = zeros(1,400);
x_0c=[0;0;0];
for i = 1:400
    phita = Vc*(Dc^i)*V_1c;
    phita_1 = Vc*(Dc^(i-1))*V_1c;
    x_k = phita*x_0c + phita_1*Bc;
    y_kc(i) = Cc*x_k;
end;
figure;
plot(y_kc);

xlabel('cicles');
ylabel('Forma controlable');
title('Forma Controlable F5');

%Observable Form

y_ko = zeros(1,400);
x_0o=[0;0;0];
for i = 1:400
    phita = Vo*(Do^i)*V_1o;
    phita_1 = Vo*(Do^(i-1))*V_1o;
    x_k = phita*x_0o + phita_1*Bo;
    y_ko(i) = Co*x_k;
end;
figure;
plot(y_ko);

xlabel('cicles');
ylabel('Forma observable');
title('Forma Observable F5');

%F6

%Controllable form

y_kc2 = zeros(1,400);
x_0c2=[10;0;0];
for i = 1:400
    phita = Vc*(Dc^i)*V_1c;
    phita_1 = Vc*(Dc^(i-1))*V_1c;
    x_k = phita*x_0c2 + phita_1*Bc;
    y_kc2(i) = Cc*x_k;
end;
figure;
plot(y_kc2);
xlabel('cicles');
ylabel('Forma controlable');
title('Forma Controlable F6');
    


%Observable Form

y_ko2 = zeros(1,400);
x_0o2=[10;0;0];
for i = 1:400
    phita = Vo*(Do^i)*V_1o;
    phita_1 = Vo*(Do^(i-1))*V_1o;
    x_k = phita*x_0o2 + phita_1*Bo;
    y_ko2(i) = Co*x_k;
end;
figure;
plot(y_ko2);
xlabel('cicles');
ylabel('Forma Observable');
title('Forma Observable F6');


%F7

[y0,t] = initial(sc, [10;0;0], 40);
[y1, t] = impulse(sc, 40);
yf = y0+y1*T;

figure;
hold on;
plot(yf)
plot(y0)
plot(y1*T)
title("Grafics F7");
axis([0 400 0 1.5]);
hold off;
legend("Suma", "Condicions Inicials", "Impuls Unitari")


%G

delta = 0.8;
ts = 1;
wn = 4/(delta*ts);
wd = wn*sqrt(1-0.8^2);
alpha1 = -2*cos(wd*T)*exp(-delta*wn*T);
alpha2 = exp(-2*delta*wn*T);
alpha3 = 0;
%G3
vector_p_z = [1 alpha1 alpha2];

vector_p_z2 = [1 alpha1 alpha2 alpha3];

a1 = den(2);
a2 = den(3);
a3 = den(4);
I=eye(3);
Kx = [(alpha3-a3) (alpha2-a2) (alpha1-a1)];

%G4

scg = ss(Ac-Bc*Kx, Bc, Cc, 0, T);
figure;
step(scg);
xlabel('segons');
ylabel('Forma controlable');
title('Forma Controlable G4');
stepinfo(scg)

%G5
Kr = inv(Cc*inv(I-Ac+Bc*Kx)*Bc);

%G6

scg2 = ss(Ac-Bc*Kx, Bc*Kr, Cc, 0, T);
figure;
step(scg2);

xlabel('segons');
ylabel('Forma controlable');
title('Forma Controlable G6');
stepinfo(scg2)

%H
%H1
CA = Cc*Ac;
A_ = [Ac(1) Ac(4) Ac(7) 0;
      Ac(2) Ac(5) Ac(8) 0;
      Ac(3) Ac(6) Ac(9) 0;
      -CA(1) -CA(2) -CA(3) 1];
CB = Cc*Bc 
B_ = [Bc(1);
    Bc(2);
    Bc(3);
    -CB];

C_ = [Cc 0];

I4 = eye(4);

syms z
matriuP = [z -1 0 0;
           0 z -1 0;
           -Ac(3) -Ac(6) (z-Ac(9)) 0;
           CA(1) CA(2) CA(3) (z-1)
           ];
determinant = det(matriuP);
coef = vpa(coeffs(determinant),5);
detMatriuP1 = [1 -Ac(9) -Ac(6) -Ac(3)];
detMatriuP = [ 1 (-Ac(9)-1) (-Ac(6)+Ac(9)) (-Ac(3)+Ac(6)) Ac(3)]
detMatriuP2 = [1 -1];
roots(detMatriuP);
A_B_ = A_*B_;
A_A_B_ = A_^2*B_;
A_A_A_B_ = A_^3*B_;
matriuCont = [B_(1) A_B_(1) A_A_B_(1) A_A_A_B_(1);
              B_(2) A_B_(2) A_A_B_(2) A_A_A_B_(2);
              B_(3) A_B_(3) A_A_B_(3) A_A_A_B_(3);
              B_(4) A_B_(4) A_A_B_(4) A_A_A_B_(4)];

% matriuPolinomi = [detMatriuP(4) detMatriuP(3) detMatriuP(2) 1;
%                   detMatriuP(3) detMatriuP(2) 1 0;
%                   detMatriuP(2) 1 0 0;
%                   1 0 0 0];
matriuPolinomi = double([coef(2) coef(3) coef(4) 1;
                  coef(3) coef(4) 1 0;
                  coef(4) 1 0 0;
                  1 0 0 0]);
Tc = matriuCont*matriuPolinomi;
Tc_1 = inv(Tc);
A_b = Tc_1*A_*Tc;


%H2
ph2 = [1 -2*cos(wd*T)*exp(-0.8*wn*T) exp(-2*0.8*wn*T) 0 0 ];

%H3
alpha1h = ph2(2);
alpha2h = ph2(3);
alpha3h = ph2(4);
alpha4h = ph2(4);

a1h = -A_b(16);
a2h = -A_b(12);
a3h = -A_b(8);
a4h = -A_b(4);

K_ = [(alpha4h-a4h) (alpha3h-a3h) (alpha2h-a2h) (alpha1h-a1h)]*Tc_1;

%H4

ssh4 = ss(A_-B_*K_, [0;0;0;1],C_,0,T);
figure;
step(ssh4);
stepinfo(ssh4);

