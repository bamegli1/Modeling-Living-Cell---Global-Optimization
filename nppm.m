%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mnguye62 MATLAB R2019b
%Inputs:
%data_xyz: (:,3) matrix of residue coordinates
%seq: bead identities  0 = nonpolar, 1 = polar, 2 = - charge, 3 = + charge
%nppm_all: harge values for polarizability of a neutral carbon reduced units, dipole moment of water, plus charge, minus
%charge, sigma, episilon for hh, episilon for hp
%k1: force constant of bead to bead bonds
%le: equilibrium length of the bonds
%temp: temperature of the NVT simulation
%zeta:  constant for BBK algorithmn
%dt: time step
%Outputs:
%U: scalar potential energy of the system 
%Fc: (:,3) matrix of all the forces in the system using more specific
%interactions than LJ potential
%interactions: are commented below
%Fr: Random force used for LD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U, Fc, Fr] = nppm(data_xyz,seq,nppm_all,k1,le,temp,m,zeta,dt)
x = data_xyz(:,1);
y = data_xyz(:,2);
z = data_xyz(:,3);

pol = nppm_all(1);
dip = nppm_all(2);
plus = nppm_all(3);
min = nppm_all(4);
sigma = nppm_all(5);
hh = nppm_all(6);
hp = nppm_all(7);

k = 1;

Fc = zeros(length(data_xyz),3);
a = zeros(length(data_xyz));

Fr = sqrt((2*temp*m*zeta)/dt)*randn(length(data_xyz),3);
count = 0;

for i = 1:size(data_xyz,1)
    for j = i+2:(size(data_xyz,1))
        
        dx = x(i,1)-x(j,1);
        dy = y(i,1)-y(j,1);
        dz = z(i,1)-z(j,1);
        
        R = sqrt(dx^2 + dy^2 + dz^2);
        R_2 = R^2;
        dxyz = [dx, dy, dz];
        Fs =  ((-12*(sigma^12)/(R_2^7))+6*(sigma^6)/(R_2^4));
        
        
        if seq(i) == 0
            if seq(j) == 0
               % induced dipole-induced dipole
                a(i,j) = 4 * hh * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = -dxyz * 4 * hh * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8));
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;    
            elseif seq(j) == 1
                % Permanent dipole - induced dipole
                a(i,j) = -(dip^2*pol)/(k^2*R^6)+4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((6*dip^2*pol/(k^2*R^8))+4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            elseif seq(j) == 2
                % Minus Charge - induced dipole
                a(i,j) = -(min^2*pol)/(2*k^2*R^4) + 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((4*min^2*pol/(2*k^2*R^6)) + 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz; 
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            else 
                %Plus Charge - induced dipole
                a(i,j) = -(plus^2*pol)/(2*k^2*R^4)+ 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((4*plus^2*pol/(2*k^2*R^6))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            end
        elseif seq(i) == 1
            if seq(j) == 0
                % Permanent dipole - induced dipole
                a(i,j) = -(dip^2*pol)/(k^2*R^6)+4 * hp * ((sigma/R^2)^6-(sigma/R^6)); 
                Fxyz = ((6*dip^2*pol/(k^2*R^8))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            elseif seq(j) == 1
                % Permanent dipole - Permanent dipole
                a(i,j) = -(2*dip^4)/((k^2)*(R^6)) + 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = -dxyz * (4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))+((-12*dip^4)/(3*(k^2)*(R^8))));
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            elseif seq(j) == 2
                % Minus Charge - Permanent dipole
                a(i,j) = -(min*dip)/(k*(R^2)) + 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((2*min*dip/(k*(R^4)))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
               
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            else 
               %Plus Charge - Permanent dipole
                a(i,j) = -(plus*dip)/(k*(R^2)) + 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((2*plus*dip/(k*(R^4)))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
              
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            end
        elseif seq(i) == 2
            if seq(j) == 0
                % Minus Charge - induced dipole
                a(i,j) = -(min^2*pol)/(2*k^2*R^4) + 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((4*min^2*pol/(2*k^2*R^6)) + 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz; 
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            elseif seq(j) == 1
                % Minus Charge - Permanent dipole 
                a(i,j) = -(min*dip)/(k*(R^2)) + 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((2*min*dip/(k*(R^4)))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
               
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            elseif seq(j) == 2
                % Minus Charge - Minus Charge
                a(i,j) = min*min/(k*R);
                Fxyz = ((min*min/(k*(R^3))))* -dxyz;
               
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            else 
                % Minus Charge - Plus Charge
                a(i,j) = min*plus/(k*R)+ 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((min*plus/(k*(R^3)))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8)))  * -dxyz;
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            end
        elseif seq(i) == 3
            if seq(j) == 0
                % Plus Charge - induced dipole
                a(i,j) = -(plus^2*pol)/(2*k^2*R^4)+ 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((4*plus^2*pol/(2*k^2*R^6))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz; 
            
            elseif seq(j) == 1
                % Plus Charge - Permanent dipole 
                a(i,j) = -(plus*dip)/(k*(R^2)) + 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((2*plus*dip/(k*(R^4)))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8))) * -dxyz;
              
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            elseif seq(j) == 2
                % Plus Charge - Minus Charge 
                a(i,j) = min*plus/(k*R)+ 4 * hp * ((sigma/R^2)^6-(sigma/R^6));
                Fxyz = ((min*plus/(k*(R^3)))+ 4 * hp * ((-12*(sigma^12)/(R^14))+6*(sigma^6)/(R^8)))  * -dxyz;
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            else 
                %Plus Charge - Plus Charge
                a(i,j) = plus*plus/(k*R);
                count = count + 1;
                Fxyz = (plus*plus/(k*(R^3))) * - dxyz;
                
                Fc(i,:) = Fc(i,:) + Fxyz;
                Fc(j,:) = Fc(j,:) - Fxyz;
            end
        end
            
    end
end

Vnb = sum(a,'all');
Vc = 0;

for i = 1:(length(data_xyz)-1)
    j = i + 1; 
    
    dx = x(i,1)-x(j,1);
    dy = y(i,1)-y(j,1);
    dz = z(i,1)-z(j,1);
    R_2 = dx^2 + dy^2 + dz^2;
    
    Vc = Vc + 0.5*k1*(sqrt(R_2)-le)^2;
    
    Fs = k1*(sqrt(R_2)-le)/sqrt(R_2);
        
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
