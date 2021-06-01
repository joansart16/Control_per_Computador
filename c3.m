clc;
clear;
close all;
%F1

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
G1pz =c2d(G, 1/10, 'zoh');

G2pz = minreal(G1pz, 0.0001);

%F5

%Impulse G2pz
impulse(G2pz);
xlim([0.0 40])
xlabel('temps');
ylabel('G2p(z)');
title('G2p(z)');

%Controllable form
Ao = [0 1 0; 0 0 1; 0.6552 -2.305 2.65];
Bo = [0; 0; 1];
Co = [0.0006173 0.0003185 0.0009183];
Do = [0];
sc = ss(Ao, Bo, Co, Do, 1/10);
[Vo, Do] = eig(Ao);
V_1o=inv(Vo);



k=1:1:400;

y_k = zeros(1,400);
x_0o=[0;0;0];
for i = 1: size(k)
    phita = Vo*(Do^i)*V_1o;
    phita_1 = Vo*(Do^(i-1))*V_1o;
    x_k = phita*x_0o + phita_1*Bo;
    y_k(i) = Co*x_k;
end;
plot(k, y_k);
    


%tm = Vo*(Do^k)*V_1o;
%CO = ctrb(sc);

%Observable Form
Ac = [0 0 0.6552; 1 0 -2.305; 0 1 2.65];
Bc = [0.0006173;0.0003185;0.0006173];
Cc = [0 0 1];
Dc = [0];
so = ss(Ac, Bc, Cc, Dc);
OB = obsv(so);
