%% methode pour initialiser les vitesses de toutes les particules en donnant
% la norme de la vitesse et le nombre de particules
function V = vitesse( v,Natome )
V=zeros(Natome,3);
% les directions sont choisies dans une boule pour assurer l'egalite de la 
% possibilite
for i=1:Natome
    x=(rand-0.5)*2;
    y=(rand-0.5)*2;
    z=(rand-0.5)*2;
    while(x^2+y^2+z^2>1)
        x=(rand-0.5)*2;
        y=(rand-0.5)*2;
        z=(rand-0.5)*2;
    end
    V(i,:)=v*[x,y,z]/sqrt(x^2+y^2+z^2);
end
end

