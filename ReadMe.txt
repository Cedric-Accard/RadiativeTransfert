Accard Cédric and Heying-Meléndrez Jhanene


 - Instructions for running the programs - 

With this read me file you'll find 5 pythons scripts all treating the scattering problem but in different manner. 
You can run all of them separatly as they're doing different things even though the last one cited will be the most complete and modular one.
For that you'll also find the requirements.txt containing the list of the packages used for this project (some in the list might be useless or unused but where already installed outside of the environment, sorry for that)

 - Scattering.py : Basic process for one photon in a spherical uniform medium. 
 You can change the radius value by changing the R value (at the very top of the program) and the density or cross section value by changing n or sigma (at the start of the while loop)

 - Multiple_Photons.py : Doing the same as the first one but now you can change the numbers of pohtons propagating inside the same medium
 The number can be changed with the "photons" parameter at the top of the program, other parameters as same as before
 
 - Unhomogenous density.py : The density model has been changed to a linearly decreasing one where the maximum value can be changed via the nc parameter in ## Density distribution ##
 Other parameter changes same as before
 
 For the program above, the origin point won't really matter but still can be changed if wanted, but to avoid overlaps (physical nonsense) the following one must take care of what value to enter even though it won't affect the process
 
 - MultipleSources.py : Add the possibility of treating multiples sources each with different density models (each can be tuned as wished inside their definition), radius and origin point
 The whole process has been implemented into a function taking the following parameter in entry and returning the following outputs : 
 
 Let's take for example one of the existing sources : 	o1, x1, y1, z1, xi1, zi1, R1 = Source(700, [1000,500,-600], 'constant')
	For the inputs :
		- Radius in meter and as an integer (here 700)
		- Coordinates of the center in the format : [x,y,z] (in meters) (here [1000,500,-600])
		- The density model choosable between 'constant', 'linear' and 'layer'
	For the outputs (7 of them) in order we have : 
		- o1 being the coordinates of the center of the sphere
		- x1 the trajectories along the x coordinate for each photon
		- y1 same for the y
		- z1 same for the z
		- xi1 the coordinates of the final position in x coordinate for telescope projection
		- zi1 same for the z
		- R1 the value of the radius previously entered
	
	You can either change the already existing sources values or create a new one (let's call it A for example) but you'll need to add the following lines for it to appear on the different plots : 
	 - oA, xA, yA, zA, xiA, ziA, RA = Source(R, [x,y,z], '.....') #to create the source, with the inputs for the functions as described as before
	 - colorsA = plt.cm._____(np.linspace(0, 1, photons)) 	#Replace the _____ by the name of an available colormap from matplotlib
	 - in #Plotting trajectories :
	      ax.plot(xA[l,],yA[l,],zA[l,], color=colorsA[l])
	 - in #Sphere Drawing : 
	 	sphereA = Sphere([oA[0], oA[1], oA[2]], RA)
		sphereA.plot_3d(ax, alpha=0.2, color='__') 		#Change the name of the color to the one you want the sphere to have
	 - in #Telescope... : 
	 	ax.scatter(xiA,ziA,  marker=".", color=colorsA)

		
		


	 
		
  