%A1
%Codi amb el que comprovam el resultat de l’enllaç tancat
C0=tf([1 10], [1]);
Gm=tf(16.34,[1 8.17 0]);
minreal(C0*Gm/(1+C0*Gm))

%A2
C=tf([1 10],1);
Gm=tf(16.34,[1 8.17 0]);
Gth=feedback(C*Gm,1);
[num, den]=tfdata(Gth, 'v');
[r,p,k]=residue(num, den)

C=tf([1 10],1);
Gm=tf(16.34,[1 8.17 0]);
Gth=feedback(C*Gm,1);
i=1;
Y=Gth*i;
tA2=0:0.001:2;
ytA2= 19.22*exp(-12.255*tA2).*cos(0.5545-3.6352*tA2);
figure;
plot(tA2,ytA2,'r')
xlabel('t(s)');
ylabel('y(t)');
title('Resposta del servomotor a un impuls unitari')


syms s
F = (10+s)*16.34/(s^2+24.51*s+163.4);
ilaplace(F)


%A3
syms s
G = (10+s)*16.34/((s^2+24.51*s+163.4)*s);
ilaplace(G)

tA3 = 0:0.001:2;
ytA3 = 1 -1.77.*exp(-12.255*tA3).*cos(1.05165-3.63*tA3);
figure;
plot(tA3, ytA3, 'r');
hold on 
xlabel('t(s)');
ylabel('y(t)');
title('Resposta a l''escalÃ³ unitari de la funciÃ³ Go(s)');
hold off

%B1
%Comprovat amb:
Gp=tf([16.34 163.4], [1 24.51 163.4 0]);
Gb = tf([0.18617],[1 0.15 0]);
[num, den]=tfdata(Gp*Gb, 'v');
[r,p,k] = residue(num, den); 

minreal( c2d(Gth*Gb, 1/20, 'zoh'), 0.0001)

Gb = tf([0.18617],[1 0.15 0]);
G_prima_z= minreal(c2d(Gth*Gb, 1/20, 'zoh'));
[z, p, k] = zpkdata(G_prima_z, 'v');

%B3
C=tf([1 -0.9925], [1 -0.25], 1/20)
rlocus(C*G_prima_z, 0:1:10000)

%B4
F = feedback(minreal(C*G_prima_z, 0.0001),1/20)

%B5
k = 0:1:50;
Y = 755*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);

figure;
T=0.05;
plot(k, Y,'r');
hold on;
xlabel('kT (segons)');
ylabel('y(kT)');
title('Esglao unitari de la funcio S(z)');
hold off;

%B6
stepinfo(Y);

%B7
k = 0:1:50;
Y = 0.1*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);

k = 0:1:50;
Y = 0.2*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);

k = 0:1:50;
Y = 0.3*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);

k = 0:1:50;
Y = 0.4*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);

k = 0:1:50;
Y = 0.5*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);

k = 0:1:50;
Y = 0.6*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);

k = 0:1:50;
Y = 0.7*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);


k = 0:1:50;
Y = 0.8*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);

k = 0:1:50;
Y = 0.9*1510*((-0.9996.^(k)*0.0011)+0.0006*0.9996.^k+2*1.3454e-3*0.5424.^k.*cos(2.3+0.1824.*k) -0.0006.*0.2496.^k);
stepinfo(Y);