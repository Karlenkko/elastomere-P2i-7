%% methode pour calculer le tenseur d'inertie en donnant la matrice des positions
function I=tenseurInertie(P)
global Natome m
I=zeros(3,3);
R=posRel(P);    % matrice des positions relatives
for k=1:Natome
    I(1,1)=I(1,1)+m*(R(k,2)^2+R(k,3)^2);    % m*(y^2+z^2)
    I(2,2)=I(2,2)+m*(R(k,1)^2+R(k,3)^2);    % m*(x^2+z^2)
    I(3,3)=I(3,3)+m*(R(k,1)^2+R(k,2)^2);    % m^(x^2+y^2)
    I(1,2)=I(1,2)-m*R(k,1)*R(k,2);  % -m*x*y
    I(2,3)=I(2,3)-m*R(k,2)*R(k,3);  % -m*y*z
    I(1,3)=I(1,3)-m*R(k,1)*R(k,3);  % -m*x*z
end
% le tenseur est symetrique
I(2,1)=I(1,2);
I(3,2)=I(2,3);
I(3,1)=I(1,3);
end
