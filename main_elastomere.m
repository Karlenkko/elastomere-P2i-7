clear all
close all
clc
%% initialisation des variables
global sigma epsilon k0 Natome m dt Niter kB

m = 1 ;	% masse d'un atome en kg
sigma = 1;  % distance ou le potentiel s'annule en m
epsilon = 1;    % porfonderur du puit du potentiel
dt = 0.005;	% pas du temps en s
Niter = 1e4;	% nombre d'iterations
kB = 1;	% cts de boltzman
T = 5;	% temperature en K
k0=30;  % raideur 
ftrac=[0 0 20];	% force de traction en N

Natome=20;   % nbr d'atome au milieu

videoflag=1;
Pini=[zeros(Natome+1,1) zeros(Natome+1,1) 1.5*(1:Natome+1)'];
defmoy=iteration(T,ftrac,videoflag,Pini)