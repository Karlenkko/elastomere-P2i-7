clear all
close all
clc
%%
%Initialisation et variables
global sigma epsilon k0 Natome m

m = 1 ;               % Masse des polymere en kg
sigma = 1;
epsilon = 1;
dt = 0.005;           % Pas d¡¯int¨¦gration en s
Niter = 1e4;          % Nombre d¡¯it¨¦rations
kB = 1;               % cts de boltzman
T = 0.01;                % temp¨¦rature en K
k0=1;
vtrac=5;
dx=vtrac*dt;

Natome=10;
P=zeros(Natome+1,3,Niter+1);
P(:,:,1)=position(Natome);

v=(3*kB*T/m)^0.5;
Vi=vitesse(v,Natome);
Vi=cancelTrans(Vi);
Vi=cancelRot(P(:,:,1),Vi);

vt=mean(sqrt(sum(Vi.^2,2)),1);
Tt=m*vt^2/3/kB;
lambda=(Tt/T)^0.5;
Vi=Vi/lambda;
%vt=mean(sqrt(sum(Vi.^2,2)),1);

sigma1 = zeros (Niter,1) ; 

V=zeros(Natome+1,3,Niter+1);
V(:,:,1)=Vi;
%%
%iteration
for i=1:Niter+1
P(Natome+1,:,i)=P(Natome+1,:,1);
V(Natome+1,:,i)=[0 0 vtrac]; 
end
P(Natome+1,3,:)=dx*[1:Niter+1];

% plot3(P(:,1,1),P(:,2,1),P(:,3,1),'-r','MarkerSize',15);
% axis([-20 20 -20 20 -30 55]);
% drawnow;
% grid on
% plot3(P(:,1,2),P(:,2,2),P(:,3,2),'-r','MarkerSize',15);
% axis([-20 20 -20 20 -30 55]);
% drawnow;


for i=1:Niter
 for k=2:Natome
P(k,:,i+1)=P(k,:,i)+dt*V(k,:,i)+dt^2/2*(forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k))/m;
V(k,:,i+1)=V(k,:,i)+dt/2*((forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k))/m+(forcetot(P(k+1,:,i+1),P(k,:,i+1),P(:,:,i+1),k)+forcetot(P(k-1,:,i+1),P(k,:,i+1),P(:,:,i+1),k))/m);
 end

P(1,:,i+1)=P(1,:,i)+dt*V(1,:,i)+dt^2/2*forcetot([0,0,0],P(1,:,i),P(:,:,i),1)/m;
V(1,:,i+1)=V(1,:,i)+dt/2*((forcetot([0,0,0],P(1,:,i),P(:,:,i),1)/m+forcetot([0,0,0],P(1,:,i+1),P(:,:,i),1))/m);

vi=mean(sqrt(sum(V(1:Natome,:,i+1).^2,2)),1);
Ti=m*vi^2/3/kB;
lambda=(Ti/T)^0.5;
V(1:Natome,:,i+1)=V(1:Natome,:,i+1)/lambda;

sigma1(i,1)=sqrt(sum(forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1).^2,2));
 
% if(mod(i,80)==0)
%  plot3(P(:,1,i),P(:,2,i),P(:,3,i),'.r','MarkerSize',15);
%  axis([-5 5 -5 5 -5 600]);
%  grid
%  drawnow;
%  end
end


for i=1:Niter
    Def(i)=(sqrt(sum(P(Natome+1,:,i).^2,2))-sqrt(sum(P(Natome+1,:,1).^2,2)))/sqrt(sum(P(Natome+1,:,1).^2,2));
end
sigma1(1,1)=sqrt(sum(forcetot(P(Natome,:,1),P(Natome+1,:,1),P(:,:,1),Natome+1).^2,2));
figure(1)
plot(Def',sigma1,'r')

L1=polyfit(Def',sigma1,1) ;
Eyoung1=L1(1)







