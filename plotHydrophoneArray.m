function plotHydrophoneArray(mic_dirs_deg, R)
% Stefan Schucker, 2023
mic_dirs_rad = mic_dirs_deg*pi/180;
Nmic = size(mic_dirs_deg,1);

hold on
% set up unit sphere information
numSphereFaces = 40;
[unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);
% radius of each sphere
spheresRadius = ones(Nmic,1)*0.05*R;
spheresX = R*cos(mic_dirs_rad(:,1)).*cos(mic_dirs_rad(:,2));
spheresY = R*sin(mic_dirs_rad(:,1)).*cos(mic_dirs_rad(:,2));
spheresZ = R*sin(mic_dirs_rad(:,2));
% for each given sphere, shift the scaled unit sphere by the
% location of the sphere and plot
for i=1:Nmic
    sphereX = spheresX(i) + unitSphereX*spheresRadius(i);
    sphereY = spheresY(i) + unitSphereY*spheresRadius(i);
    sphereZ = spheresZ(i) + unitSphereZ*spheresRadius(i);
    h = surf(sphereX, sphereY, sphereZ);
    set(h, 'FaceColor','white', 'FaceLighting', 'gouraud')
    
    text(1.1*spheresX(1),1.1*spheresY(1),1.1*spheresZ(1), 'Top','color','r','FontSize',24)
    text(1.1*spheresX(2),1.1*spheresY(2),1.1*spheresZ(2), 'A','color','r','FontSize',24)  
    text(1.1*spheresX(3),1.1*spheresY(3),1.1*spheresZ(3), 'B','color','r','FontSize',24)  
    text(1.1*spheresX(4),1.1*spheresY(4),1.1*spheresZ(4), 'C','color','r','FontSize',24)  
end
sphereX = unitSphereX*0.98*R;
sphereY = unitSphereY*0.98*R;
sphereZ = unitSphereZ*0.98*R;
h = surf(sphereX, sphereY, sphereZ);
set(h, 'FaceColor','green', 'FaceAlpha', 0.3)
light('Position',[0 0 1],'Style','infinite');
material dull
hold off
axis equal
grid on
view(-30,15)
set(gca,'visible','off')
