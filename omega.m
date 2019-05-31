%% calcul du vitesse angulaire en donnant le moment cinetique et le tenseur d'inertie 
function w=omega(L,I)
% L=I^omega
w=I\L;
end
