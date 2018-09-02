%--------------------------------------------------------------------------
% final version of dmft code for bachelorthesis
%
%   'Berechnung der Spektralfunktion in DMFT-Behandlung des Hubbard-Modells
%    mittels exakter Diagonalisierung'
%
% by Mirco Hellwig
%
%--------------------------------------------------------------------------
% Main dmft loop:
%       1. calculating local and bath green function using self-energy
%       2. fitting anderson model
%       3. solving anderson model using exact diagonalization or lanczos
%       4. gaining new self-energy
%       5. repeat until difference between old and new self-energy becomes
%          lower than the convergence condition
%--------------------------------------------------------------------------

% general parameters
a=1;                                            % gitter constant 
t=0.25;                                         % hopping in tight-binding/hubbard model
eps_k=tightbinding(a,t);                        % dispersion in tight-binding 
U=0:0.1:5;                                      % Hubbard correlation U, if a vector is given dmft loops over all entries
beta=1000;                                      % defining temperatur T with beta=1/tT
n_max=1000;                                     % number of matsubara frequencies
n_s=4;                                          % number of bath sites
iw_n=1i.*(2.*(0:n_max-1)+1).*pi./beta;          % matsubara frequencies i*w_n 
mixing_parameter=1;                             % determine mixing of new and old self-energy, 1 for using pure new, smaller for higher portion of old

%--------------------------------------------------------------------------
% parameters for the ED-solver
usealgorithm='ed';                              % Algorithm to use, 'ed' for exact diagonalization or 'lancz' for lanczos algortihm
if strcmp(usealgorithm,'lancz')                 % additional parameters for lanczos:
    lancziter=5000;                             % number of iterations           
    lanczstates=30;                             % number of used states
end
omegamin=-8;                                    % intervall and number of
omegamax=8;                                     % points on the real
Nomega=8096;                                    % frequency axis
broadening=1e-2;                                % broadening used for ED solver and calculating spectral function
                                               
%--------------------------------------------------------------------------
% starting loop over all given Hubbard-U's
for ii=1:length(U)
    
    %----------------------------------------------------------------------
    % Initialisation, creating directory, printing important parameters to console
    fprintf('\n----------------Init----------------\n')
    selfen=0.*iw_n;                             % starting point for self-energy Sigma(iw_n) 
    mu=U(ii)/2;                                 % chemical potential
    step=0;                                     % iteration step in dmft loop
    diff=1;                                     % difference between new and old self-energy, also used as convergence criterium
    directory=strcat('/results/beta_',num2str(beta),'_ns_',num2str(n_s),'_U_',num2str(U(ii))); % directory for results 
    
        % Catching 'Directory already exists' to prevent unintended mixing 
        % of old and new results in the same directory (2nd calculation for
        % same parameters). For some reason I don't quite understand, it's 
        % not possible to catch warnings in Matlab, so I had to convert it
        % into an error.
        s=warning('error','MATLAB:MKDIR:DirectoryExists'); %#ok<CTPCT>
        try
            mkdir(directory)
        catch
            inp_str=input(['The directory already exists. \n',...
                    'Do want to continue and delete old results? y/n \n'],'s');
            if strcmp(inp_str,'y')
                rmdir(directory,'s')
                mkdir(directory)
            else
                continue
            end
        end
        warning(s)
    
    mkdir(strcat(directory,'/cmdout'))          % directory for saving terminal output of ED-solver 

    fprintf('U=%f \nbeta=%f \nnmax=%u \nns=%u \n',[U(ii) beta n_max n_s])
    fprintf(strcat('algorithm=',usealgorithm,'\n'))
    
%--------------------------------------------------------------------------
% starting main loop
    while diff>1e-10                            % convergence condition 
        step=step+1;
        fprintf('\n---------------STEP-%u---------------\n',step)
        
        % calculate local green function, bath green function and fit anderson model
        xmin=andersonfit(mu,n_s,iw_n,selfen,eps_k);      
        [v_l,eps_l]=para(xmin);             % parameters of AIM
        
        % writing input file for ED-solver, taking used algorithm into
        % account
        if strcmp(usealgorithm,'lancz')
            writeconfig(U(ii),n_s,v_l,eps_l,beta,omegamin,omegamax,Nomega,...
                        broadening, n_max, directory, lancziter, lanczstates)
        elseif strcmp(usealgorithm,'ed')
            writeconfig(U(ii),n_s,v_l,eps_l,beta,omegamin,omegamax,Nomega,...
                        broadening, n_max, directory)
        else
            error('No such algorithm known!!')       
        end
        
        % calling ED-solver
        selfen_new=ED(directory,step);
        
        % number of electrons per atom calculated using impurity spectral
        % function as control number, should always be 1.00000
        % small deviation <0.1 seems to be ok
        electron_per_atom=numinatom(directory, step);
        fprintf('\n Electrons in atom: %5f\n',electron_per_atom)
        
        % calculate difference between new and old self-energy
        diff=sum(abs(selfen-selfen_new).^2)./numel(iw_n);
        fprintf('\n Difference: %5e\n',diff)
        
        % overwrite self-energy with the new one using defined mixing
        % parameter
        selfen=mixing_parameter.*selfen_new+(1-mixing_parameter).*selfen;
    end
    
    %----------------------------------------------------------------------
    % main loop finished
    fprintf('\n--------------CONVERGED--------------\n')
    delete_unnecessary(directory)           % clear unimportant output files of ED-solver
                                            % there are a lot!!
end

