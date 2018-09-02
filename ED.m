function selfen=ED(directory,step)
% calling the Fortran ED-routines, using the input file in given directory and
% writing all output in it
%
%   Args:
%       direcotry: directory for input and output
%       step:      iteration step of dmft loop
%
%   Returns:
%       selfen:    new self-energy on Matsubara frequencies

fprintf('\n->Calling fortran ED-routine<-\n')

working_directory=cd(directory);

% directory of the ED-solver
ED_directory='~/Dropbox/Bachelorarbeit/Programme/ED_ast11';

%-------------------------------------------------------------------------
% calling Fortran ED routines
[~,cmdout] = system(strcat(ED_directory,'/eigen'));
dlmwrite(strcat('cmdout/eigen_',num2str(step)),cmdout,'');

[~,cmdout] = system(strcat(ED_directory,'/selfen'));
dlmwrite(strcat('cmdout/selfen_',num2str(step)),cmdout,'');

[~,cmdout] = system(strcat(ED_directory,'/bundle2plot selfen_matsub_bundle.dat'));
dlmwrite(strcat('cmdout/bundle2plot_',num2str(step)),cmdout,'');

%--------------------------------------------------------------------------
% obtaining new self-energy
A=dlmread('selfen_matsub_1_1.dat');
selfen=(A(:,3)+1i.*A(:,4))';

cd(working_directory)

