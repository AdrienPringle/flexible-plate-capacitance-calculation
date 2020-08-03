import femm
import math
import matplotlib
import matplotlib.pyplot as plt

#######################################################################
# Large parts of this code are adapted from David Meeker's example code
# available at femm.info/wiki/pyFEMM
#######################################################################

# data setup

#number of trials from 0 - 1 flexpercent (total number will be *2 - 1 to account for both sides)
TRIALCOUNT = 5

# simulation setup

#number of nodes on top plate
NODECOUNT = 20
#plate radius
RADIUS = 18
#dielectric width
DI_WIDTH = 0.05
#base distance between plates
DISTANCE = 0.5
#space between dielectric and top plate
SPACE = DISTANCE - DI_WIDTH
# r distance between nodes
dr = RADIUS/NODECOUNT

def getPlateZ(flexPercent, r):
	if (r == 0) : return SPACE * flexPercent
	jVar = math.pi * r / RADIUS #input to the J_0 function
	return SPACE * flexPercent * math.sin(jVar) / jVar

def analyseCapacitor(flexPercent):

	# Create a new electrostatics problem
	femm.newdocument(1)

	# Draw the geometry
	femm.ei_probdef('centimeters','axi',10**(-8),30)
	femm.ei_drawrectangle(0,-SPACE,RADIUS,-DISTANCE) # bottom plate and dielectric
	# femm.ei_drawline(0,0,RADIUS,0) # top plate
	for i in range(0,NODECOUNT):
		femm.ei_drawline(i*dr, getPlateZ(flexPercent, i*dr), (i+1)*dr, getPlateZ(flexPercent, (i+1)*dr))

	femm.ei_drawarc(0,-RADIUS * 5,0,RADIUS * 5,180,2.5) # arc boundary
	femm.ei_drawline(0,-RADIUS * 5,0,RADIUS * 5) # line boundary

	# Create and assign a "fixed voltage" boundary condition to curve
	femm.ei_addboundprop('Fixed Voltage',0.5,0,0,0,0)
	femm.ei_selectarcsegment(0,-RADIUS*5)
	femm.ei_setarcsegmentprop(2.5,'Fixed Voltage',0,0,'<none>')
	femm.ei_clearselected()

	# Add and assign the block labels for the air and dielectric regions
	femm.ei_addmaterial('air',1,1,0)
	femm.ei_addmaterial('dielectric',3,3,0)

	femm.ei_addblocklabel(10,10)
	femm.ei_selectlabel(10,10)
	femm.ei_setblockprop('air',0,1,0)
	femm.ei_clearselected()

	femm.ei_addblocklabel(RADIUS/2,-DISTANCE + DI_WIDTH/2)	
	femm.ei_selectlabel(RADIUS/2,-DISTANCE + DI_WIDTH/2)
	femm.ei_setblockprop('dielectric',0,1,0)
	femm.ei_clearselected()

	# Add a "Conductor Property" for each of the plates
	femm.ei_addconductorprop('anode',1,0,1)
	femm.ei_addconductorprop('cathode',0,0,1)

	# Assign the anode properties to all sides of the first strip
	# femm.ei_selectsegment(RADIUS/2,0)
	for i in range(0, NODECOUNT):
		femm.ei_selectsegment((i+0.5)*dr, getPlateZ(flexPercent, (i+0.5)*dr))

	femm.ei_setsegmentprop('<None>',0.25,0,0,0,'anode')
	femm.ei_clearselected()

	# Assign the cathode properties to all sides of the fourth strip
	femm.ei_selectsegment(RADIUS/2,-DISTANCE)
	femm.ei_setsegmentprop('<None>',0.25,0,0,0,'cathode')
	femm.ei_clearselected()

	femm.ei_zoomnatural()

	# Save the geometry to disk so we can analyze it
	femm.ei_saveas('capacitors/capacitorflex'+ str(flexPercent) + '.fee')

	femm.ei_analyze()
	femm.ei_loadsolution()

	femm.eo_selectblock(10,10)
	energy = femm.eo_blockintegral(0)[0]


	return energy

# Start up and connect to FEMM
femm.openfemm()

energies = [2*analyseCapacitor(i/TRIALCOUNT) for i in range(-TRIALCOUNT, TRIALCOUNT+1)]

# close femm after all simulations are finished 
femm.closefemm()

x = [i/TRIALCOUNT for i in range(-TRIALCOUNT, TRIALCOUNT+1)]

fig, ax = plt.subplots()
ax.plot(x,energies, 'o-')

#set plot style
ax.set(
	xlabel='Percent of Total Displacement (-1 is closer together, 1 is farther apart)',
	ylabel='Capacitance (F)',
	title='Impact of Drum Membrane Flex on Capacitance'
)
plt.grid()

#show plot
plt.show()
