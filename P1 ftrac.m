clc
close all
clearvars
%%
%Initialisation et variables
global sigma epsilon k0 Natome m

m = 1 ;               % Masse des polymere en kg
sigma = 1;
epsilon = 1;
dt = 0.005;           % Pas d¡¯int¨¦gration en s
Niter = 1e4;          % Nombre d¡¯it¨¦rations
kB = 1;               % cts de boltzman
T = 500;                % temp¨¦rature en K
k0=400;
vtrac=5;
dx=vtrac*dt;
ftrac=[0,0,1000];

Natome=20;
P=zeros(Natome+1,3,Niter+1);
P(:,:,1)=position(Natome);

v=(3*kB*T/m)^0.5;
Vi=vitesse(v,Natome);
Vi=cancelTrans(Vi);
Vi=cancelRot(P(:,:,1),Vi);
E=zeros(Niter,1);

vt=(sum(Vi(1,:).^2,2))^0.5;
Tt=m*vt^2/3/kB;
lambda=(Tt/T)^0.5;
Vi=Vi*lambda;

sigma1 = zeros (Niter,1) ;   
%%
%iteration
P(:,:,2)=P(:,:,1)+dt*Vi;
% for i=1:Niter+1
%     P(Natome+1,3,i)=dx*i;
% end

plot3(P(:,1,1),P(:,2,1),P(:,3,1),'-r','MarkerSize',15);
axis([-20 20 -20 20 -30 55]);
drawnow;
grid on
% plot3(P(:,1,2),P(:,2,2),P(:,3,2),'-r','MarkerSize',15);
% axis([-20 20 -20 20 -30 55]);
% drawnow;

%%
for i=2:Niter
%     P(Natome+1,:,i+1)=2*P(Natome+1,:,i)-P(Natome+1,:,i-1)+dt^2/m*(force(P(Natome,:,i),P(Natome+1,:,i))+ftrac);   
%  for k=2:Natome
%      P(k,:,i+1)=2*P(k,:,i)-P(k,:,i-1)+dt^2/m*(force(P(k+1,:,i),P(k,:,i))+force(P(k-1,:,i),P(k,:,i)));
%  end
 
 
     P(Natome+1,:,i+1)=2*P(Natome+1,:,i)-P(Natome+1,:,i-1)+dt^2/m*(forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1)+ftrac);   
 for k=2:Natome
     P(k,:,i+1)=2*P(k,:,i)-P(k,:,i-1)+dt^2/m*(forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k));
 end
 
 
 fi=force(P(Natome,:,i),P(Natome+1,:,i));
 
%  P(Natome+1,1,i+1)=2*P(Natome+1,1,i)-P(Natome+1,1,i-1)+dt^2/m*fi(1);
%  P(Natome+1,2,i+1)=2*P(Natome+1,2,i)-P(Natome+1,2,i-1)+dt^2/m*fi(2);
%  
% P(1,:,i+1)=2*P(1,:,i)-P(1,:,i-1)+dt^2/m*(force([0,0,0],P(1,:,i))+force(P(2,:,i),P(1,:,i)));
  P(1,:,i+1)=2*P(1,:,i)-P(1,:,i-1)+dt^2/m*(forcetot([0,0,0],P(1,:,i),P(:,:,i),1)+forcetot(P(2,:,i),P(1,:,i),P(:,:,i),1));
  
 sigma1(i,1)=(sum(force(P(Natome,:,i),P(Natome+1,:,i)).^2,2))^0.5;
 if(mod(i,80)==0)
 plot3(P(:,1,i),P(:,2,i),P(:,3,i),'.r','MarkerSize',15);
 axis([-5 5 -5 5 -5 100]);
 grid
 drawnow;
 end
end
for i=1:Niter
    Def(i,1)=(sum((P(Natome+1,:,i)-P(Natome+1,:,1)).^2,2))^0.5;
end
sigma1(1,1)=(sum(force(P(Natome,:,1),P(Natome+1,:,1)).^2,2))^0.5;
figure(1)
plot(Def,sigma1,'r')
mean(Def)
temp=polyfit(Def,sigma1,1);
sigma2=temp(1)*Def+temp(2);
figure(3)
plot(Def,sigma2,'r')
% V=zeros(Natome+1,3,Niter+1);
% V(:,:,1)=Vi;
% for i=2:Niter
% V(:,:,i)=(P(:,:,i+1)-P(:,:,i-1))/2/dt;
% E(i)=sum(sum(V(:,:,i).^2*m*0.5,2),1);
% end
% E(1)=sum(sum(Vi.^2*m*0.5,2),1);
% figure(2)
% plot(E);






