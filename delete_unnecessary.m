function delete_unnecessary(directory)
% delete all unnecessary output files of the ED-solver after dmft-convergence
%
%   Args:
%       directory:  directory with said output files

for f=dir(directory).'
    if (~f.isdir && ~any(strcmp(f.name,{'dos.dat','selfen_1_1.dat','selfen_matsub_1_1.dat','hubbard1.cfg'})))
                                        % ^ the only important ones ^
        delete(fullfile(directory, f.name))
    end
end

