%% methode pour calculer les positions relatives par rapport au centre de masse 
% d'un ensemble des particules en donnant la matrice des positions
function R=posRel(P)
global Natome
posG=centerDeMasse(P);  % centre de masse 
R=zeros(Natome,3);  % matrice des positions relatives
for k=1:Natome
    R(k,:)=P(k,:)-posG(1,:);    % positions relatives
end
end

