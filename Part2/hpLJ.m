%mnguye62 MATLAB R2019b
function [U, Fc, Fr] = hpLJ(data_xyz,seq,hp_all,sigma,k,le,temp,m,zeta,dt)
x = data_xyz(:,1);
y = data_xyz(:,2);
z = data_xyz(:,3);

hh = hp_all(1);
hp = hp_all(2);
pp = hp_all(3);

Fc = zeros(length(data_xyz),3);
a = zeros(length(data_xyz));

Fr = sqrt((2*temp*m*zeta)/dt)*randn(length(data_xyz),3);

for i = 1:length(data_xyz)
    for j = i+2:(length(data_xyz))

        dx = x(i,1)-x(j,1);
        dy = y(i,1)-y(j,1);
        dz = z(i,1)-z(j,1);
        R_2 = dx^2 + dy^2 + dz^2;
        Fs =  ((-12*(sigma^12)/(R_2^7))+6*(sigma^6)/(R_2^4));
        
        
            if seq(i) == 'H'
                if seq(j) == 'P'
                    
                    a(i,j) = 4 * hp * ((sigma/R_2)^6-(sigma/R_2^3));
        
                    Fx = -dx *4 * hp * Fs;
                    Fy = -dy *4 * hp * Fs;
                    Fz = -dz *4 * hp * Fs;
                    
                    Fc(i,1) = Fc(i,1) + Fx;
                    Fc(i,2) = Fc(i,2) + Fy;
                    Fc(i,3) = Fc(i,3) + Fz;
        
                    Fc(j,1) = Fc(j,1) - Fx;
                    Fc(j,2) = Fc(j,2) - Fy;
                    Fc(j,3) = Fc(j,3) - Fz;
                    
                else

                    a(i,j) = 4*hh*((sigma/R_2)^6-(sigma/R_2^3));
        
                    Fx = -dx *4 * hh * Fs;
                    Fy = -dy *4 * hh * Fs;
                    Fz = -dz *4 * hh * Fs;
                    
                    Fc(i,1) = Fc(i,1) + Fx;
                    Fc(i,2) = Fc(i,2) + Fy;
                    Fc(i,3) = Fc(i,3) + Fz;

                    Fc(j,1) = Fc(j,1) - Fx;
                    Fc(j,2) = Fc(j,2) - Fy;
                    Fc(j,3) = Fc(j,3) - Fz;
                end
            else
                if seq(j) == 'P'

                    a(i,j) = 4*pp*((sigma/R_2)^6-(sigma/R_2^3));
        
                    Fx = -dx *4 * pp * Fs;
                    Fy = -dy *4 * pp * Fs;
                    Fz = -dz *4 * pp * Fs;
                    
                    Fc(i,1) = Fc(i,1) + Fx;
                    Fc(i,2) = Fc(i,2) + Fy;
                    Fc(i,3) = Fc(i,3) + Fz;

                    Fc(j,1) = Fc(j,1) - Fx;
                    Fc(j,2) = Fc(j,2) - Fy;
                    Fc(j,3) = Fc(j,3) - Fz;
                else
                    
                    a(i,j) = 4*hp*((sigma/R_2)^6-(sigma/R_2^3));
        
                    Fx = -dx *4 * hp * Fs;
                    Fy = -dy *4 * hp * Fs;
                    Fz = -dz *4 * hp * Fs;
                    
                    Fc(i,1) = Fc(i,1) + Fx;
                    Fc(i,2) = Fc(i,2) + Fy;
                    Fc(i,3) = Fc(i,3) + Fz;

                    Fc(j,1) = Fc(j,1) - Fx;
                    Fc(j,2) = Fc(j,2) - Fy;
                    Fc(j,3) = Fc(j,3) - Fz;
                end
            end
    end
end

Vnb = sum(a(:));
Vc = 0;

for i = 1:(length(data_xyz)-1)
    j = i + 1; 
    
    dx = x(i,1)-x(j,1);
    dy = y(i,1)-y(j,1);
    dz = z(i,1)-z(j,1);
    R_2 = dx^2 + dy^2 + dz^2;
    
    Vc = Vc + 0.5*k*(sqrt(R_2)-le)^2;
    
    Fs = k*(sqrt(R_2)-le)/sqrt(R_2);
        
    Fx = -dx * Fs;
    Fy = -dy * Fs;
    Fz = -dz * Fs;
                    
    Fc(i,1) = Fc(i,1) + Fx;
    Fc(i,2) = Fc(i,2) + Fy;
    Fc(i,3) = Fc(i,3) + Fz;

    Fc(j,1) = Fc(j,1) - Fx;
    Fc(j,2) = Fc(j,2) - Fy;
    Fc(j,3) = Fc(j,3) - Fz; 
end

U = Vc + Vnb;
