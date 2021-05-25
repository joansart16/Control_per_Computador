K=0.5;
J=0.009;
R=3.4;
r=0.0125;
L=0.47;
d=0.15;
Km=150;
a=0.5*(K^2/(J*R));
Gm=tf(K/(J*R), [1 K^2/(J*R) 0]);
C0 = tf([Km Km*a], [1]);
Gb = tf([(5/7)*9.81*(r/L)],[1 d 0]);
Gth=feedback(C0*Gm,1);
Gth=minreal(Gth,0.0001);
G = minreal(Gth*Gb, 0.01);

%F2