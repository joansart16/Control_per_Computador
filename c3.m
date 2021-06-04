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

%F5

%Impulse G2pz

impulse(G2pz);
xlim([0.0 40])
xlabel('temps');
ylabel('G2p(z)');
title('G2p(z)');

[num, den] = tfdata(G2pz, 'v');

%Controllable form
Ac = [0 1 0; 0 0 1; -den(4) -den(3) -den(2)];
Bc = [0; 0; 1];
Cc = [num(4) num(3) num(2)];
Dc = 0;
sc = ss(Ac, Bc, Cc, Dc, T);
[Vc, Dc] = eig(Ac);
V_1c=inv(Vc);





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
    


%tm = Vo*(Do^k)*V_1o;
%CO = ctrb(sc);

%Observable Form
Ao = [0 0 -den(4); 1 0 -den(3); 0 1 -den(2)];
Bo = [num(4);num(3);num(2)];
Co = [0 0 1];
Do = 0;
so = ss(Ao, Bo, Co, Do, T);

[Vo, Do] = eig(Ao);
V_1o=inv(Vo);




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
    


%Observable Form

y_ko2 = zeros(1,400);
x_0o2=[10;0;0];
for i = 1:400
    phita = Vo*(Do^i)*V_1o;
    phita_1 = Vo*(Do^(i-1))*V_1o;
    x_k = phita*x_0o2 + phita_1*Bo;
    y_ko2(i) = Co*x_k;
end;
plot(y_ko2);


%F7

[y0,t] = initial(sc, [10;0;0], 40);
[y1, t] = impulse(sc, 40);
yf = y0+y1*T;
figure;
plot(yf);



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

Kx = [(alpha3-a3) (alpha2-a2) (alpha1-a1)];



%G5
Kr = inv(Cc*inv([1 0 0; 0 1 0; 0 0 1] -Ac + Bc*Kx)*Bc);





