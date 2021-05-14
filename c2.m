K=0.5;
J=0.009;
R=3.4;
r=0.0125;
L=0.47;
d=0.15;
Gm=tf(K/(J*R), [1 K^2/(J*R) 0]);
c=tf([150, (150*(K^2)/(J*R))], 1);
Gth=feedback(c*Gm,1);
Gth=minreal(Gth,0.0001);
Gb = tf([(5/7)*9.81*(r/L)],[1 d 0]);
G1 = c2d(Gth*Gb, 1/5);
G2=zpk(minreal(G1,0.001));
%figure(1); clf; hold on; zgrid; rlocus(G2);
syms x
wn = solve((1.9704-2*cos(x*sqrt(1-0.7^2)*(1/5))*exp(-0.7*x*(1/5)))/0.0036753== (exp(-2*0.7*x*(1/5))-0.9704)/(0.003675*0.9981)  ,x);
 
K1=double((1.9704-2*cos(wn*sqrt(1-0.7^2)*(1/5))*exp(-0.7*wn*(1/5)))/0.0036753);
C1=tf([K1],[1]);
S1=minreal(feedback(K1*G1, 1),0.0001);

%C5
step(S1);
xlabel('temps(s)');
ylabel('S1(z)');
title('Resposta del sistema S1(z)');
stepinfo(S1)

%C6a
syms x
wn2 = double(solve( (-2*cos(x*sqrt(1-0.7^2)/5)*exp(-0.7*x*1/5)+1)/0.0036753 == exp(-2*0.7*x*1/5)/(0.9981*0.0036753) ,x));
 
K2 = double(exp(-2*0.7*wn2*1/5)/(0.9981*0.0036753))

%C6b
Td = 0.9704/0.148;
K2p = K2/(1+Td*5)


