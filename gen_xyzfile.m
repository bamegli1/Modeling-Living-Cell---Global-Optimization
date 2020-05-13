%Read in a list of coordinates, formated into 3 columns, such that
%x1 y1 z1
%x2 y2 z2
%etc. Each configuration of Np particles is immediately followed by another
%configuration at a new timestep. See, e.g. shortconfig.dat
%
%load the file. Here the matlab variable it reads in is named tot_conf.
%tot_conf is a numerical array, with 3 columns and Np*Number_of_configs rows.
%Np is the number of particles per configuration.

function[]=gen_xyzfile(tot_conf, Np, varargin)

for i = 1:2:length(varargin)
    name = varargin{i};
    this = varargin{i+1};
    switch lower(name)
        case 'skip'
            skip = this;
        case {'file', 'filename', 'fname'}
            fname = this;
    end
end
if ~exist('skip', 'var')
    skip = 1;
end
if ~exist('fname', 'var')
    fname = 'crdsVMD';
end

[L, col]=size(tot_conf);
Nc=L/Np;

Str='H';
iter='Iteration';
fid=fopen(strjoin({fname, '.xyz'}, ''),'w');
t=1;
for i=1:skip:Nc
    fprintf(fid,'%d\n',Np);
    fprintf(fid,'%s\n',iter);
    for j=1:1:Np
        fprintf(fid,'%s %2.3f %2.3f %2.3f\n',Str, tot_conf(t,1),tot_conf(t,2),tot_conf(t,3));
        t=t+1;
    end
    t = t+(skip-1)*Np;
end

fclose(fid);
