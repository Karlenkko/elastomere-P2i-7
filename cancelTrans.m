%% methode pour annuler la vitesse de translation d'un ensemble de particules
% apres l'initialisation aleatoire des vitesse en donnant la matrice des vitesses
function newvit = cancelTrans (V)
Vtrans=mean(V,1);   % moyenne de la vitesse en m
newvit=V;
% annulation de la vitesse de translation selon 3 axes
for i=1:size(V,2)
    newvit(:,i)=V(:,i)-Vtrans(1,i); 
end
end