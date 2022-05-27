% Test time length 
tfinal = 540; %mins 
time = 1:tfinal;
I = size(time);

% Relax stage
I(1:540) = 60;

% Internal Resistance (Ohms)
Rin = 0.00154;
% Tab Resistance (Ohms)
Rtab = 9.74e-5;
% Heat Flux Total
Qgen=size(time);
for t=time 
    Qgen(t)=(I(t)^2)*Rin;
end

% Plot Current vs Time
figure
plot(time,I)
title('Current Test Schedule','interpreter','latex','Fontsize',14)
xlabel('Time (mins)','interpreter','latex','Fontsize',14)
ylabel('Current (A)','interpreter','latex','Fontsize',14)
% Plot Qin vs Time
figure
plot(time,Qgen)
title('Heat in','interpreter','latex','Fontsize',14)
xlabel('Time (mins)','interpreter','latex','Fontsize',14)
ylabel('Heat flux (W)','interpreter','latex','Fontsize',14)