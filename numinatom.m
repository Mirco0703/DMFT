function electron_per_atom=numinatom(directory, step)
% reads the number of electrons per atom that the ED-solver calculates
% using the spectral density
%
%   Args:
%       directory:          directory with output files of ED-solver
%
%       Returns:
%       electron_per_atom:  number of electrons per atom

fid=fopen(strcat(directory,'/cmdout/selfen_', num2str(step)));
electron_per_atom=0;
% reads each line of output file (->variable length) and compare it to
% 'Number of particles...', number is given by characters 58 till end in this line
while electron_per_atom==0
    str=fgetl(fid);
    k=strfind(str,'Number of particles in the atom from spectral density:');
    if k
        electron_per_atom=str2double(str(58:end));
    end
end
    