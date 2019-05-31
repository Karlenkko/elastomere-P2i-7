%% calcul du vitesse angulaire en donnant le moment cinetique et le tenseur d'inertie 
function w=omega(L,I)
% L=I^omega
% probleme de 0/0=NaN
for i=1:size(I,1)
    for j=1:size(I,2)
        if I(i,j)==0
            I(i,j)=eps;
        end
    end
end
w=I\L;
end
