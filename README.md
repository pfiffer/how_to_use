Instructions to use the Direct Optimization Method (PSO) on SPE-10 model

Title: New Efficient Method for Injection Well Location Optimization using Fast Marching Method

Author: Reza Yusefzadeh

Email: reza_yusef@yahoo.com

Files included:

1- PSO.m: Comprised of main configuration sections of the PSO algorithm and parameters setting.
2- NPV.m: Takes well positions from PSO.m and assings it to ECLDATA variable and sends it to Eclipse.m
3- Eclipse.m: Takes well positions from NPV.m and writes them to "sched.m" on 5th line after the "WELSPECS". Then calls ECLIPSE software using dos('$multirun.bat') command.
   Finally, extracts FOPT, FWPT, and FWIT from "SPE10.data" and send them back to NPV.m to calculate the Net Present Value
4- Perm-het.txt: Conatains the heterogenous permeability values.
5- SPE10.data: Model characterization information.
6- sched.dat: Contains well locations
7- perm.inc: Permeability distribution file.
8- Perm-het.txt: Permability distribution to be loaded from MATLAB.
9- FMM_2.m: Fast Marching Method to calculate Time of Flight and Production wells drainage boundary.

Permeability map is loaded by MATLAB at line 16of Main.m
You should provide filename containing permeability values in the argument of load() function. This file should only contain numeric values.
Number of colmns in permeability file is entered in line 23. The default value is 6.

Production well locations are specified in line 42 to be used to calculate drainage boundaries.

Then you should set the number of Maximum Iterations and number of Population denoted by MaxIt and nPop respectively.

When you run the prgram, you should enter the number of decision which is the double the number of wells to be optimized.
Then enter number of production wells previously present in the model. Default is 6 for SPE-10 model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Start of the Fast Marching Method Function

Time of Flight and drainage area of each production wells are assigned to TOF and Index variables resepetively in line 44.
Afterwards, regions where indices from wells are different are determined using Index variable to determine the drainage boundary.
This region is taken as the search domain of the optimization algorithm.

TOF and Index is calculated using "FMM_2" function. This function takes the production well positions and permeability as the input data.
Next, you should specify the number and dimension of grid blocks of the model and the other properties of the fluid and rock to calculate the diffusivity constant.
Index of each region is calculated and assigned to FPI (Frozen Points Index) variable and returned.

End of the Fast Marching Method Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To change the number of wells to be optimized (Default is 2), first of all, comment or uncomment the lines pretaining to that wells in the "sched.dat" file.

You can change the Oil price, drilling, water injection and production costs in "NPV.m" file.

Results are saved every iteration in line 226.


Output files:
Three output files are provided for test. Black squares on the figures indicate the results from PSO method.

1- "Well location for SPE10 (6 prods and 2 injs).fig"

2- "Well location for SPE10 (6 prods and 3 injs).fig"

3- "Well location for SPE10 (6 prods and 4 injs).fig"


Software Requirements:
		You need to have Eclipse 2010.1 and MATLAB R2013a or higher versions to run these files.

Note: If you are using MATLAB R2015a, you need to change the 48th line of the Eclipse.m to:
""""""""""""""""	AA = cell2mat(A.textdata(4));	""""""""""""""""""""""
(Only the text, not quotations and space)
