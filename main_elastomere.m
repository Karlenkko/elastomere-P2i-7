clear all
close all
clc
%% initialisation des variables
global sigma epsilon k0 Natome m dt Niter kB

m = 1 ;	% masse d'un atome en kg
sigma = 1;  % distance ou le potentiel s'annule en m
epsilon = 1;    % profonderur du puit du potentiel
dt = 0.005;	% pas du temps en s
Niter = 3e4;	% nombre d'iterations
kB = 1;	% cts de boltzman
T = 5;	% temperature en K
k0=30;  % raideur 
ftrac=[0 0 0];	% force de traction en N

Natome=20;   % nbr d'atome au milieu

videoflag=1;
Pini=[zeros(Natome+1,1) zeros(Natome+1,1) 1.5*(1:Natome+1)'];
% defmoy=iteration(T,ftrac,videoflag,Pini)

%% demo libre
defmoy=iteration(T,ftrac,videoflag,Pini);

%%
Def=zeros(3,3);
videoflag=0;
figure(1)
for i=1:size(Def,1)
    T=i*2.5;
    for j=1:size(Def,2)
        ftrac=[0 0 15+5*j];
        Def(i,j)=iteration(T,ftrac,videoflag,Pini);     
    end
    plot(Def(i,:),20*(1:size(Def,2)));
    axis([-5 50 0 70]);
    hold on
end