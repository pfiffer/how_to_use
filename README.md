# how_to_use
Direct.m: Main loops of the algorithm which sends well locations to the NPV.m \n
NPV.m: Gets well locations and sends them to Eclipse.m \n
Eclipse.m: Contains Eclipse() function which takes well locations ans writes them to the sched.dat file then calls ECLIPSE software to simulate the model with new well configuration \n
Perm-het.txt: Containig heterogeneous permeability values.\n
HTCM.data: containig characteristics of the synthetic models (Included for test). \n

Note: The original file is arranged for heterogeneous cases. If you want to change it to a homogenous one, you should make the following changes to some file:
1- Uncomment lines 158 to 160 in HTCM.data file.
2- As the NPV map will be symmetric, it will be enough to calculate only a quarter of the NPV map and finally replicate it to other parts
	So, you should change the upper bounds of the "for" loops to "13".

3- Also, uncomment lines 35 and 35, then comment lines 45, 47-49 and uncomment line 46.


Software Requirements:
		You need to have Eclipse 2010.1 and MATLAB R2013a or higher versions to run these files.

Note: If you are using MATLAB R2015a, you need to change the 47th line of the Eclipse.m to:
""""""""""""""""	AA = cell2mat(A.textdata(4));	""""""""""""""""""""""
(Only the text, not quotations and space)
