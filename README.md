# ospats
Julia package for optimal spatial stratification

The package ospats contains 3 Julia script files:
1. “main”: the script for reading the grid data, setting process parameters and calling the function ospatsmr() or function allocatemr()  
2. “ospatsmr”: containing only the function ospatsmr().
3. “allocatemr”: containing only the function allocatemr(). This function should be used if coarse-gridding is necessary, 
i.e. if the grid has too many points to enable processing of the pairwise distance matrix.

The functions ospatsmr() and allocatemr() both lead to a stratification by the Ospats method. 
The difference is that opspatsmr() applies the method to the entire grid, while allocatemr() applies the method to a systematic sample 
from the grid, followed by allocation of the remaining grid points to the strata as constructed from the sample. 

The requested function is called by setting the parameter "in" (from sampling INterval) in the main script. 
If in=1 then ospatsmr() is called, if “in” is set to a value larger than 1, then allocatemr() is called. 
For instance, if in=3, then every third grid point gets included in the sample.

For further explanation see the user manual at https://github.com//jjdegruijter/ospats/Ospats-manual.

