%% methode qui donne la position du centre de masse des particules en donnant
% la matrice des positions ou les particules ont la meme masse
function posG=centerDeMasse(P)
global Natome m
posG=zeros(1,3);
% addition des position selon 3 axes
for k=1:Natome
    posG(1,1)=posG(1,1)+(P(k,1)*m); % selon x
    posG(1,2)=posG(1,2)+(P(k,2)*m); % selon y
    posG(1,3)=posG(1,3)+(P(k,3)*m); % selon z
end
posG=posG/Natome/m;  
end