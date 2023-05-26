#!/bin/tcsh -f

set Lengths = (0 1 2 3)

foreach Length ($Lengths)
	echo "Length = $Length"

	./MonteCarlo $Length
	
end
