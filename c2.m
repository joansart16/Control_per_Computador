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

%C7
C2p = K2p*tf([(1+(Td/(1/5))) -Td*5], [1 0], 1/5);
S2 = minreal(feedback(C2p*G1, 1), 0.0001)

%C8

step(S2);
xlabel('temps');
ylabel('S2(z)');
title('Resposta del sistema S2(z)');
stepinfo(S2)


%D2
syms a
g = double(solve(4/(a*2) == pi/(1.2*sqrt(1-a^2)), a));
wn3 = 4/(g*2);
Grs=tf([wn3^2],[1 2*wn3*g wn3^2]);
r=roots([1 2*wn3*g wn3^2]);
t2=-exp(r(1)/5)-exp(r(2)/5);
t3=exp(r(1)/5)*exp(r(2)/5);
Grz = c2d(Grs ,1/5, 'matched')


%D5
C3=minreal((1/G2)*Grz/(1-Grz))

%D6
S3=minreal(feedback(C3*G1, 1), 0.0001)

%D7
step(S3);
xlabel('temps(s)');
ylabel('S3(z)');
title('Resposta del sistema S3(z)');
stepinfo(S3)

%D8



%E
G1e = c2d(Gth*Gb, 1);
G2e=minreal(zpk(G1e), 0.00001);


%E1



%E3

K4=1/(0.08862*0.9628 +0.08862);
b0 = K4*0.08862*0.9628;

%E4
C4=tf([K4 -K4*0.8607], [1 b0], 1);

%E5
S4 = zpk(minreal(feedback(C4*G1e, 1), 0.0001));

%E6
step(S4);
ylim([-0.1 1.2])
xlabel('temps');
ylabel('S2(z)');
title('Resposta del sistema S4(z)');
stepinfo(S4)

%E7

