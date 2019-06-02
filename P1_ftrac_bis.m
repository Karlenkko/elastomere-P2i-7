clear all
close all
clc
%% initialisation des variables
global sigma epsilon k0 Natome m

m = 1 ;	% masse d'un atome en kg
sigma = 1;  % distance ou le potentiel s'annule en m
epsilon = 1;    % porfonderur du puit du potentiel
dt = 0.005;	% pas du temps en s
Niter = 1e4/2;	% nombre d'iterations
kB = 1;	% cts de boltzman
T = 1;	% temperature en K
k0=30;  % raideur 
ftrac=[0 0 20];	% force de traction en N

Natome=20;   % nbr d'atome au milieu

%% initialisation des positions
P=zeros(Natome+1,3,Niter+1);    % positions des particules
P(:,:,1)=[zeros(Natome+1,1) zeros(Natome+1,1) 1.5*(1:Natome+1)'];
sigma1 = zeros (Niter,1) ; 
%% calcul de vitesse aleatoire
v=(3*kB*T/m)^0.5;
Vi=vitesse(v,Natome+1);
Vi=cancelTrans(Vi);
Vi=cancelRot(P(:,:,1),Vi);
%% rescaling des vitesses
vt=mean(sqrt(sum(Vi.^2,2)),1);
Tt=m*vt^2/3/kB;
lambda=(Tt/T)^0.5;
Vi=Vi/lambda;
vt=mean(sqrt(sum(Vi.^2,2)),1);
V=zeros(Natome+1,3,Niter+1);
V(:,:,1)=Vi;
%% iteration du dernier atome
% for i=1:Niter+1
%     P(Natome+1,:,i)=P(Natome+1,:,1)+dx*(i-1);
%     V(Natome+1,:,i)=vtrac; 
% end
figure(1)
%% iteration pour tous
% schema de Verlet avec un pas
% P(n+1)=P(n)+dt*V(n)+dt^2/2*F(n)
% V(n+1)=V(n)+dt/2*(F(n+1)+F(n))
for i=1:Niter
%% calcul de P
    for k=2:Natome
        fextcurrent=forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k);
        P(k,:,i+1)=P(k,:,i)+dt*V(k,:,i)+dt^2/2*fextcurrent/m;
%         V(k,:,i+1)=V(k,:,i)+dt*fextcurrent/m;        
    end
    fext1current=forcetot([0,0,0],P(1,:,i),P(:,:,i),1)+forcetot(P(2,:,i),P(1,:,i),P(:,:,i),1);
    P(1,:,i+1)=P(1,:,i)+dt*V(1,:,i)+dt^2/2*fext1current/m;

%     V(1,:,i+1)=V(1,:,i)+dt/2*fext1current/m;
    
    fextdercurrent=forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1)+ftrac;
    P(Natome+1,:,i+1)=P(Natome+1,:,i)+dt*V(Natome+1,:,i)+dt^2/2*fextdercurrent/m;
%     V(Natome+1,:,i+1)=V(Natome+1,:,i)+dt*fextdercurrent/m;

%% calcul de V
    for k=2:Natome
        fextcurrent=forcetot(P(k+1,:,i),P(k,:,i),P(:,:,i),k)+forcetot(P(k-1,:,i),P(k,:,i),P(:,:,i),k);
        fextnext=forcetot(P(k+1,:,i+1),P(k,:,i+1),P(:,:,i+1),k)+forcetot(P(k-1,:,i+1),P(k,:,i+1),P(:,:,i+1),k);
        V(k,:,i+1)=V(k,:,i)+dt/2*(fextcurrent/m+fextnext/m);        
    end
    fext1next=forcetot([0,0,0],P(1,:,i+1),P(:,:,i+1),1)+forcetot(P(2,:,i+1),P(1,:,i+1),P(:,:,i+1),1);
    V(1,:,i+1)=V(1,:,i)+dt/2*(fext1current/m+fext1next/m);    
    
    fextdernext=forcetot(P(Natome,:,i+1),P(Natome+1,:,i+1),P(:,:,i+1),Natome+1)+ftrac;
    V(Natome+1,:,i+1)=V(Natome+1,:,i)+dt/2*(fextdercurrent/m+fextdernext/m);
 
%%    
%     thermosta
    vi=mean(sqrt(sum(V(1:Natome,:,i+1).^2,2)),1);
    Ti=m*vi^2/3/kB;
    lambda=sqrt(Ti/T);
    V(1:Natome,:,i+1)=V(1:Natome,:,i+1)/lambda;

    sigma1(i,1)=sqrt(sum(forcetot(P(Natome,:,i),P(Natome+1,:,i),P(:,:,i),Natome+1).^2,2));
    
%     video
	if(mod(i,10)==0)
    	Ptemp=[zeros(1,3);P(:,:,i)];  
        plot3(Ptemp(:,1),Ptemp(:,2),Ptemp(:,3),'.-r','MarkerSize',25);
%         view(0,0);
        axis([-10 10 -10 10 -40 60]);
        grid
        drawnow;
	end
end

for i=1:Niter
%     Def(i)=(sqrt(sum(P(Natome+1,:,i).^2,2))-sqrt(sum(P(Natome+1,:,1).^2,2)))/sqrt(sum(P(Natome+1,:,1).^2,2));
%     Def(i)=sqrt(sum(P(Natome+1,:,i).^2,2)
    Def(i)=sqrt(sum(P(Natome+1,:,i).^2,2))-sqrt(sum(P(Natome+1,:,1).^2,2));
end
sigma1(1,1)=sqrt(sum(forcetot(P(Natome,:,1),P(Natome+1,:,1),P(:,:,1),Natome+1).^2,2));
figure(2)
plot(Def',sigma1,'r')
defmoy=mean(Def)
L1=polyfit(Def',sigma1,1) ;
Eyoung1=L1(1)

for i=1:Niter+1
vtt(i)=mean(sqrt(sum(V(1:Natome,:,i).^2,2)),1);
Ttt(i)=m*vtt(i)^2/3/kB;
end
figure(3)
plot(Ttt);
figure(4)
plot(vtt);




