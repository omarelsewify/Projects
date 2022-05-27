        % Hole
    for i = 49:71
        for j = 65:89
            k_hole = 0.133;
            rho_hole = 867;
            cp_hole = 1945;
            alpha2 = k_hole/(rho_hole*cp_hole);
            Fo2 = (alpha2*dt)/(gs^2);            % Fourier Number Coolant
            phi(i,j) = phi_old(i,j)+Fo2*(phi_old(i,j+1)-4*phi_old(i,j) +phi_old(i,j-1) +phi_old(i+1,j) +phi_old(i-1,j));
        end
    end
    
    term2 = (k_hole/k);
    
    % Hole Left
    for i = 48
        for j = 64:88
            term2 = (k_hole/k);
            phi(i,j) = Fo*(phi_old(i,j+1) +2*phi_old(i-1,j) +phi_old(i,j-1) +2*term2*phi_old(i+1,j) +((1/Fo) -4 -(2*term2))*phi_old(i,j));
        end
    end
    
    % Hole Right
    for i = 72
        for j = 64:88
            phi(i,j) = Fo*(phi_old(i,j+1) +2*phi_old(i+1,j) +phi_old(i,j-1) +2*term2*phi_old(i-1,j) +((1/Fo) -4 -(2*term2))*phi_old(i,j));
        end
    end
    
    % Hole Bottom
    for i = 48:72
        for j = 64
            phi(i,j) = Fo*(phi_old(i-1,j) +2*phi_old(i,j-1) +phi_old(i+1,j) +2*term2*phi_old(i,j+1) +(1/Fo -4 -2*term2)*phi_old(i,j));
        end
    end
    
    % Hole Top
    for i = 48:72
        for j = 88
            phi(i,j) = Fo*(phi_old(i-1,j) +2*phi_old(i,j+1) +phi_old(i+1,j) +2*term2*phi_old(i,j-1) +((1/Fo) -4 -(2*term2))*phi_old(i,j));
        end
    end