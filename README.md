# how_to_use

Instructions to use the PSO-FMM method on the Synthetic models


Title: New Efficient Method for Injection Well Location Optimization using Fast Marching Method

Author: Reza Yusefzadeh

Email: reza_yusef@yahoo.com

Files included:

1- Direct.m: Comprised of main loops connecting to NPV.m and the plotting sections
2- NPV.m: Takes well positions from Direct.m and assings it to ECLDATA variable and sends it to Eclipse.m
3- Eclipse.m: Takes well positions from NPV.m and writes them to "sched.m" on 5th line after the "WELSPECS". Then calls ECLIPSE software using dos('$multirun.bat') command.
   Finally, extracts FOPT, FWPT, and FWIT from "HTCM.data" and send them back to NPV.m to calculate the Net Present Value
4- Perm-het.txt: Conatains the heterogenous permeability values.
5- sched.dat: Contains well locations.
6- Perm-het.txt: Permeability distribution file to be loaded from MATLAB.
12- FMM_2: Function to calculate Time of Flight and drainage boundaries.


Enter the number of grid blocks in x- and y-direction in line 15. Default is 25 for the both directions.
To convert the permeability values to a matrix representing the model, two "for" loops are used. Specify the number of columns in the permeability file in line 21. Default is 25.


****** Enter the location of production wells in line 49 to calculate the Time of Flight and drainage boundaries. ********

Note: The original file is arranged for heterogeneous cases. If you want to change it to a homogenous one, you should make the following changes to some file:
1- Uncomment lines 158 to 160 in HTCM.data file.
2- Comment lines 17 to 29 and uncomment lines 33 to 35 in Main.m file

3- Also, uncomment lines 35 and 35, then comment lines 45, 47-49 and uncomment line 46.


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



Output files:
Three output files are provided to test the drainage boundaries. Black squares on the figures indicate the location of production wells. 

1- "Drainage boundary of 4 prod wells heterogeneous case.fig"

2- "Drainage boundary of 4 prod wells homogeneous case.fig"

Software Requirements:
		You need to have Eclipse 2010.1 and MATLAB R2013a or higher versions to run these files.

Note: If you are using MATLAB R2015a, you need to change the 47th line of the Eclipse.m to:
""""""""""""""""	AA = cell2mat(A.textdata(4));	""""""""""""""""""""""
(Only the text, not quotations and space)
