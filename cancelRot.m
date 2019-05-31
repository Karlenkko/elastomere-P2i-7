%% methode pour annuler la vitesse de rotation apres l'initialisation aleatoire 
% des vitesseen donnant la matrice des positions et la matrice des vitesses
function vSR=cancelRot(P,V)
global Natome
vSR=V;
R=posRel(P);    % positions relatives
I=tenseurInertie(P);    % tenseur d'inertie
L=momentInertie(P,V);   % moment cinetique
w=omega(L,I);   % vitesse angulaire
% annulation des vitesses de rotation selon 3 axes
% -omega^R
for k=1:Natome
    vSR(k,1)=vSR(k,1)+(w(2,1)*R(k,3)-w(3,1)*R(k,2));    % composante x        
    vSR(k,2)=vSR(k,2)+(w(3,1)*R(k,1)-w(1,1)*R(k,3));    % composante y
    vSR(k,3)=vSR(k,3)+(w(1,1)*R(k,2)-w(2,1)*R(k,1));    % composante z
end
end