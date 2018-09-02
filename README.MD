--------------------------------------------------------------------------

 final version of dmft code for bachelorthesis

   'Berechnung der Spektralfunktion in DMFT-Behandlung des Hubbard-Modells
    mittels exakter Diagonalisierung'

 by Mirco Hellwig

--------------------------------------------------------------------------

 Main dmft loop:
 	1. calculating local and bath green function using self-energy
    2. fitting anderson model
    3. solving anderson model using exact diagonalization or lanczos
    4. gaining new self-energy
    5. repeat until difference between old and new self-energy becomes lower than the covergence condition

--------------------------------------------------------------------------

Installation:
	1. working ED-solver of the Arbeitsgruppe Wehling, ITP Uni Bremen
		-> change ED_directory in ED.m accordingly
	2. working Matlab(R)
		I personally have used 2015b, but it should be pretty indepentant of the exact version
	3. that should be it, simply run dmft.m

--------------------------------------------------------------------------

*  Parameters are all hardcoded at the beginning of dmft.m, I have not made the effort of an extra input 	file
*  Mott transition happens at to small U (approx. a factor 2), but the reason could not be found
	
