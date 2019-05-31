%% methode pour calculer le vecteur du moment cinetique en donnant 
% la matrice des vitesses et la matrice des positions 
function L=momentInertie(P,V)
global Natome m
L=zeros(3,1);
R=posRel(P);    % matrice des positions relatives 
% L=somme(vecteur position^quantite deu mouvement)
for k=1:Natome
    L(1,1)=L(1,1)-m*(R(k,2)*V(k,3)-R(k,3)*V(k,2));  % composante x
    L(2,1)=L(2,1)-m*(R(k,3)*V(k,1)-R(k,1)*V(k,3));  % composante y
    L(3,1)=L(3,1)-m*(R(k,1)*V(k,2)-R(k,2)*V(k,1));  % composante z
end
end