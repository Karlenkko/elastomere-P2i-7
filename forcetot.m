%% methode pour calculer la force de FENE+LJ soumise d'une molecule, 
% en donnant la position d'elle-meme, la molecule voisine (pour calculer
% FENE),la matrice de position de l'ensemble, et le numero de la molecule
function [F] = forcetot(A,B,P,n)
%%¨force A sur B (FENE)
global sigma epsilon k0
r=((A(1)-B(1))^2+(A(2)-B(2))^2+(A(3)-B(3))^2)^0.5;  % module de distance
e=[(A(1)-B(1))/r,(A(2)-B(2))/r,(A(3)-B(3))/r];  % vecteur de force
f=k0*(r-1.5*sigma); % moudule de force de FENE
F=f*e;  % force sous forme vectorielle

%% force de LJ
rc=1;
coeff=1;
for i=1:size(P,1)
    if(i~=n)
        rtemp=((P(i,1)-B(1))^2+(P(i,2)-B(2))^2+(P(i,3)-B(3))^2)^0.5;    % module de distance
        if(rtemp<rc)
        etemp=[(P(i,1)-B(1))/rtemp,(P(i,2)-B(2))/rtemp,(P(i,3)-B(3))/rtemp];    % vecteur de force
        F=F-(0.5*coeff*24*epsilon/sigma*(2*(sigma/rtemp)^13-(sigma/rtemp)^7))*etemp;    % force sous forme vectorielle
        end
    end
end
% force de LJ de la molecule sur l'autre dans le  point fixe
r0=((B(1))^2+(B(2))^2+(B(3))^2)^0.5;
if(r0<rc)
e0=[(0-B(1))/rtemp,(0-B(2))/rtemp,(0-B(3))/rtemp];
F=F-(0.5*coeff*24*epsilon/sigma*(2*(sigma/r0)^13-(sigma/r0)^7))*e0;
end
end