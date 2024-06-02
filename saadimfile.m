m1=5;m2=1;l=0.2;g=10;k=0.5;b1=0.25;b2=0.25;
N=50;R=11.6;J=0.006;fm=0.0009;km=0.954;
K=km/(N*R);B=fm+(km^2/R);
D22=m1+m2+N^2*J;D12=m2*l;D11=m2*l^2+N^2*J;
D=[D11 D12;D12 D22];
M=inv(D)
M11=M(1,1);
M12=M(1,2);
M21=M(1,2);

M22=M(2,2);



syms x1 x2 x3 x4
[X1,X2,X3,X4]=solve(x2==0,-M11*((b1+N^2*B)*x2-m2*g*l*sin(x1))-M12*((b1+N^2*B)*x4+k*x3-m2*l*x2^2*sin(x1))==0,x4==0,-M12*((b1+N^2*B)*x2-m2*g*l*sin(x1))-M22*((b1+N^2*B)*x4+k*x3-m2*l*x2^2*sin(x1))==0)

J=jacobian([x2;-M11*((b1+N^2*B)*x2-m2*g*l*sin(x1))-M12*((b1+N^2*B)*x4+k*x3-m2*l*x2^2*sin(x1));x4;-M12*((b1+N^2*B)*x2-m2*g*l*sin(x1))-M22*((b1+N^2*B)*x4+k*x3-m2*l*x2^2*sin(x1))])

Alinear=subs(J,[x1 x2 x3 x4],[X1 X2 X3 X4])

A1=[0 1 0 0;M11*m2*g*l -M11*(b1+N^2*B) -M12*k -M12*(b1+N^2*B);0 0 0 1;M12*m2*g*l -M12*(b1+N^2*B) -M22*k -M22*(b1+N^2*B)];
eigA1=eig(A1)
B1=[0 0;M11*N^2*K M12*N^2*K;0 0;M12*N^2*K M22*N^2*K];
C1=[1 0 0 0;0 0 1 0];
D1=[0 0;0 0];
sys1=ss(A1,B1,C1,D1)
sys=tf(sys1)

t=0:0.1:10;
u=[0*t;0*t];
[ty,yt,X]=lsim(sys1,u,t,[2 2 2 2]);
figure(1)
plot(X(:,1),X(:,2))
figure(2)
plot(X(:,3),X(:,4))

figure(3)
bode(sys1)

figure(4)
nyquist(sys1)

figure(5)
step(sys1)

figure(6)
pzmap(sys1)