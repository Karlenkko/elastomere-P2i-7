%% methode d'initialisation des positions des atomes entourees par une boule 
% unitaire en donnant le nombre de particules
function P = position( Natome )
P=zeros(Natome,3);
for i=1:Natome
    x=(rand-0.5)*2;
    y=(rand-0.5)*2;
    z=(rand-0.5)*2;
%     contrainte de boule
    while(x^2+y^2+z^2>1)
        x=(rand-0.5)*2;
        y=(rand-0.5)*2;
        z=(rand-0.5)*2;
    end
    P(i,:)=[x,y,z]/((x^2+y^2+z^2)^0.5);
end
end

