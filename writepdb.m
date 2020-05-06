%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% writepdb
% Melissa Mai
% April 2017
%
% Creates a PDB file for a linear chain of bonded atoms for visualization
% in PyMOL. Only useful for two different beads.
%
% INPUTS
% config: (npart*nconfig)x3 matrix of coordinates for each atom at each
%           time step
% numid: npartx1 vector of IDs for each bead. 
%           1 = P (polar), 0 = H (hydrophobic)
% step: stepsize. eg, set to 1 to visualize every time step or set to 100
%           to visualize every 100th configuration.
% fileID: filename (without the extension)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writepdb(config, numid, step, fileID)
    % Define file name
    if ~exist('fileID', 'var')
        fileID = input('Input name of structure: ', 's');
    end
    filename = strjoin({fileID, '.pdb'}, '');
    pdbfile = fopen(filename, 'w');
    
    % Determine number of particles
    npart = size(numid,1);
    
    % Total number of configurations
    nconfig = size(config,1)/npart;
    xc = reshape(config(:, 1), npart, nconfig);
    yc = reshape(config(:, 2), npart, nconfig);
    zc = reshape(config(:, 3), npart, nconfig);
    
    % Filter by step
    xc = xc(:, 1:step:nconfig);
    yc = yc(:, 1:step:nconfig);
    zc = zc(:, 1:step:nconfig);
    
    nstruct = size(xc,2);
    
    for i = 1:nstruct
        fprintf(pdbfile, 'MODEL\n');
        for j = 1:npart
            % Choose O and C for visualization in PyMOL
            if numid(j) == 1
                aname = 'O';
                resname = 'POL';
            else
                aname = 'C';
                resname = 'HYD';
            end
            
            occ = 1.00; bfac = 0.00;
            
            % PDB format as of 2017?
            fprintf(pdbfile, '%s  %5d  %-3s%4s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n', 'ATOM', j, aname, resname, j, xc(j,i), yc(j,i), zc(j,i), occ, bfac);
        end
        for j = 1:npart-1
            fprintf(pdbfile, '%s%5d%5d\n', 'CONECT', j, j+1);
        end
        % End config
        fprintf(pdbfile, 'ENDMDL\n');
    end
    % End file
    fprintf(pdbfile, 'END');
    fclose(pdbfile);