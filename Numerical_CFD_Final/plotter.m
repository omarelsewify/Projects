function plotter(V,xc,yc)
gam = 1.4;

close all

figure
surf(xc,yc,V(:,:,1))
title_text=('Density vs Position');
title(title_text,'interpreter','latex','Fontsize',14);
xlabel('x (m)','interpreter','latex','Fontsize',14)
ylabel('y (m)','interpreter','latex','Fontsize',14);
zlabel('Density (kg/m3)','interpreter','latex','Fontsize',14)
colorbar
view(0,90)


figure
surf(xc,yc,V(:,:,2))
title_text=('VelocityX vs Position');
title(title_text,'interpreter','latex','Fontsize',14);
xlabel('x (m)','interpreter','latex','Fontsize',14)
ylabel('y (m)','interpreter','latex','Fontsize',14);
zlabel('VelocityX (m/s)','interpreter','latex','Fontsize',14)
colorbar
view(0,90)

figure
surf(xc,yc,V(:,:,3))
title_text=('VelocityY vs Position');
title(title_text,'interpreter','latex','Fontsize',14);
xlabel('x (m)','interpreter','latex','Fontsize',14)
ylabel('y (m)','interpreter','latex','Fontsize',14);
zlabel('VelocityY (m/s)','interpreter','latex','Fontsize',14)
colorbar
view(0,90)

figure
surf(xc,yc,V(:,:,4))
title_text=('Pressure vs Position');
title(title_text,'interpreter','latex','Fontsize',14);
xlabel('x (m)','interpreter','latex','Fontsize',14)
ylabel('y (m)','interpreter','latex','Fontsize',14);
zlabel('Pressure (Pa)','interpreter','latex','Fontsize',14)
colorbar
view(0,90)


figure
for i = 1:60
    for j = 1:20
        c(i,j) = sqrt(gam*V(i,j,4)./V(i,j,1));
        M(i,j) = V(i,j,2)/c(i,j);
    end
end
surf(xc,yc,M(:,:))
title_text=('Mach Number vs Position');
title(title_text,'interpreter','latex','Fontsize',14);
xlabel('x (m)','interpreter','latex','Fontsize',14)
ylabel('y (m)','interpreter','latex','Fontsize',14);
zlabel('Mach','interpreter','latex','Fontsize',14)
colorbar
view(0,90)

end