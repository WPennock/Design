# Bhalunka Sedimentation Tank Design
Created by William Pennock on 7 August 2020 for the Bhalunka 1 L/s plant in collaboration with AguaClara Reach and Gram Vikas. The calculations were made with reference to:

- The [textbook](https://aguaclara.github.io/Textbook/Sedimentation/Sed_Intro.html)
- The design Mathcad worksheet
- The [CEE 4520 slides](https://github.com/monroews/CEE4520Lectures)
- An incomplete version of [AIDE](https://github.com/AguaClara/aguaclara/blob/master/aguaclara/design/sed_tank.py)

## Imports
```python
from aguaclara.play import *
from aguaclara.core import pipes
```

## General Parameters
```python
# Plant
FlowPlant = 1*(u.L/u.s) # Design flow rate
VolStorage = 6000*u.L # Nominal storage of water tower
TimeRunPlant = (VolStorage/FlowPlant).to(u.hr) # Single run time to fill water tower.
print('The time for a single plant run is ',np.floor(TimeRunPlant.magnitude),' hour and ',np.round((TimeRunPlant.magnitude%1)*60,1),' minutes.',sep='')

# Physical
Temperature = 15*u.degC # Average low temperature for Bhubaneswar
nu = pc.viscosity_kinematic_water(Temperature) # Kinematic viscosity of water
rho = pc.density_water(Temperature) #  Density of water
```

## Sedimentation Tank Design
 This is a typical masonry sedimentation tank, but with a few key differences:
1. To eliminate the "lost triangle" under the last plate settler and keep the average velocity through the plate settlers closer to the upflow velocity, an inner wall will be built within the sedimentation tank at the location of the bottom of the last plate settler. This will direct flow up at a higher rate and will serve as the weir for the floc hopper, with the floc hopper occupying the space between the inner wall (weir) and the outer wall. Originally, a 60° sloped wall was proposed, but Santiago at APP thought this would be less practical.
2. Having one sedimentation tank removes the necessity of the four-channel system. In its place, there will be an influent pipe from the pipe flocculator to the inlet manifold and a single pipe from outlet manifold to disinfection.
3. To make it possible to drain poorly flocculated water before it reaches the sedimentation tank, it is important to have a tee with a ball valve on the inlet pipe.
4. An upflow velocity of 0.85 mm/s (as opposed to the standard 1 mm/s) was selected as a conservative design in the absence of knowledge of the composition of the raw water. High organic content water could have low density flocs that will have a lower settling velocity that requires a floc blanket with a correspondingly lower upflow velocity.

**Questions:**
- How is the sed tank water level controlled in the plantitas (i.e., without the control weir)?
- How should the drain for the floc hopper be sized?

### General Parameters
Primarily using typical AguaClara conventions, as noted in the code. The upflow velocity is lower to account for the likely presence of natural organic matter in flocs.
```python
VelUpSed = 0.85*u.mm/u.s # Upflow velocity of sedimentation tank
AngleSlope = 50*u.deg # Angle of sloped walls at bottom of sedimentation tank
WidthSed = 1*u.m # Width of sedimentation tank (this is standard)
LengthSedMin = FlowPlant/(VelUpSed*WidthSed) # Length of sedimentation tank
LengthSedMin.to(u.m)
LengthSed = 1.2*u.m
VelUpSed = FlowPlant/(LengthSed*WidthSed)

DepthFlocBlanket = 1*u.m # ? This sets the height of the sedimentation tank.
## In addition to providing sufficient collisions, this also needs to be big enough to make flow uniform into the plate settlers.
```

### Floc Characteristics
This is open for debate. For design purposes, I am still running with assumptions of kaolin flocs with perhaps some adsorption of DOM. We didn't analyze the composition of the water, so we only know it can get up to around 10 NTU, and this turbidity was presumably both organic and inorganic.

```python
ShearFlocMax = 0.5*u.Pa # The maximum allowable shear for flocs, based on the textbook derivation from Garland
DensityCompactFloc = 2*u.kg/u.m**3 # The density of flocs after they have settled and compacted in the floc hopper. This is an unsubstantiated number just to get an estimate.
```

### Inlet Diffuser Design
```python
# Inputs
WidthDiffuser = 1/8*u.inch # Mold fabricated with 1/8" plate steel
HeightDiffuser = 15*u.cm # The diffusers are 15 cm long pieces of pipe
HeadLossDiffuser = 1*u.cm # Standard AguaClara design assumption
DiamNomDiffuser = 1*u.inch # Standard nominal diameter of diffusers
PiStretchDiffuser = 1.2 # Assumed 20% stretch when creating flare
# SpaceDiffuser = 6*u.cm # B_Diffuser (was 6 cm in plantita)

# Calculations

## Maximum Diffuser Velocity
# By minor loss equation
VelMaxDiffuser = (np.sqrt(2*u.g_0*HeadLossDiffuser)).to(u.mm/u.s)
# Calculated based on minor loss equation with K = 1.
print('The maximum velocity of the sedimentation tank diffusers by head loss is ',VelMaxDiffuser)
# By estimated maximum shear on flocs
PiPlaneJet = 0.0124 # Taken from textbook
VelMaxDiffuserShear = ((ShearFlocMax/rho)**(1/2)*(VelUpSed*WidthSed/(nu*PiPlaneJet))**(1/4)).to(u.mm/u.s) # Eq. 434
print('The maximum velocity of the sedimentation tank diffusers by fluid shear is ',VelMaxDiffuserShear,'.',sep='')
# Diffuser Velocity
VelDiffuser = (VelUpSed*WidthSed/WidthDiffuser).to(u.mm/u.s) # This jet velocity is lower than the head loss maximum and the shear maximum, so this appears to be adequate. WidthJetReversed = WidthSed*VelUpSed/VelDiffuser
print('The sedimentation diffuser velocity is ',VelDiffuser,'.',sep='')
WidthJetReversed
EDRDiffuser = (PiPlaneJet*(VelDiffuser**3/WidthJetReversed)).to(u.mW/u.kg)
HeadLossDiffuser = (VelDiffuser**2/(2*u.g_0)).to(u.cm)
print('The head loss across the diffusers is ',HeadLossDiffuser,'.',sep='')

## Diffuser Flare Dimensions
LocPipeDiffuser = (np.abs(np.array(pipes.pipedb['NDinch'])-DiamNomDiffuser.magnitude)).argmin() # Index for diffuser in pipe_database.csv
AreaPVCDiffuser = (np.pi/4)*((pipes.pipedb.iloc[LocPipeDiffuser,1]*u.inch)**2-(pipes.pipedb.iloc[LocPipeDiffuser,6]*u.inch)**2) # Cross sectional area of 1" PVC pipe (Schedule 40)

ThickStretchDiffuser = pipes.pipedb.iloc[LocPipeDiffuser,5]*u.inch/PiStretchDiffuser # Stretched wall thickness of diffuser
LengthOutDiffuser = (AreaPVCDiffuser/2/ThickStretchDiffuser - WidthDiffuser).to(u.cm)# Diffuser Outer Length, B_Diffuser
LengthInDiffuser = (LengthOutDiffuser - (2*ThickStretchDiffuser)).to(u.cm) # Diffuser inner length, S_Diffuser
print('The outer length of the diffuser is ',LengthOutDiffuser, ', and the inner length is ',LengthInDiffuser,'.',sep='')

## Number of Diffusers
NumDiffuser = np.floor((LengthSed/LengthOutDiffuser).to(u.dimensionless)) # Assumes that diffusers will be touching each other.
print('The number of diffusers is ',NumDiffuser.magnitude,'.',sep='')
```
The inlet diffusers will be 15 cm long segments of 1" pipe molded such that they have a flare that is about 5.3 cm long (outer length). These will be placed close together to approximate a continuous line jet, requiring 22 diffusers. They have a sufficiently low velocity to not exceed the design head loss or the maximum shear on the flocs. They have a head loss of about 3.5 cm.

### Inlet Manifold Design
```python
PiFlowManifold = 0.8
VelMaxManifold = (VelDiffuser*np.sqrt(2*(1-PiFlowManifold**2)/(PiFlowManifold**2+1))).to(u.mm/u.s) # Eq. 435
DiamMinManifold = pc.diam_circle(FlowPlant/VelMaxManifold).to(u.inch)
DiamMinManifold
LocPipeManifold = (np.array(pipes.pipedb['ID_SCH40'])>DiamMinManifold.magnitude).argmax() # Index for manifold in pipe_database.csv
DiamNomManifold = pipes.pipedb.iloc[LocPipeManifold,0]
print('The nominal diameter of the inlet manifold is ',DiamNomManifold,' inches.',sep='')

HeightReverser2Diffuser = 3*u.cm
```
The calculated inlet manifold diameter is 3.5" with a length equal to the length of the sedimentation tank bay (1.2 m). If that diameter is not available, 4" pipe will work. The jet reverser nominal diameter is assumed to be 3" with a length of 1.2 m.

### Plate Settler Design
```python
AnglePlate = 60*u.deg
SpacePlate = 2.5*u.cm
VelCapture = 0.12*u.mm/u.s
ThickPlate = 2*u.mm # This is an assumption of single wall thickness

LengthMinPlate = ((SpacePlateSed*(VelUpSed/VelCaptureSed-1)+ThickPlateSed*(VelUpSed/VelCaptureSed))/(np.sin(AnglePlateSed)*np.cos(AnglePlateSed))).to(u.cm)
LengthPlate = np.ceil(LengthMinPlate)
print('The length of the plate settlers is ',LengthPlate,sep='')

HeightPlate = np.sin(AnglePlate)*LengthPlate
print('The height of the plates in the sedimentation tank is ',np.round(HeightPlate,1),'.',sep='')
OverlapHorizontalPlate = np.cos(AnglePlate)*LengthPlate
print('The horizontal overlap of the plates is ',OverlapHorizontalPlate,'.',sep='')

SpaceHorizontalPlate = ((ThickPlate + SpacePlate)/np.sin(AnglePlate)).to(u.cm)
NumPlate = np.floor((LengthSed/SpaceHorizontalPlate).to(u.dimensionless))
print('The number of sedimentation plates required is ',NumPlate,'.',sep='')
```
The plate settlers will be spaced 2.5 cm apart at a 60° angle, resulting in 38 plates 0.38 m in length. Within the sedimentation tank, these will occupy a height of about 33 cm. The plate modules will be assembled with 1/2" PVC with 3/4" spacers and will be supported with 1.5" PVC. The plates overlap a horizontal distance of 19 cm, which means that extending the length of the sedimentation tank by 20 cm will accommodate the overlap without the need of angled walls (the floc weir will help eliminate the "lost triangle".)

### Floc Hopper Design
The floc hopper will be built within the walls of the sedimentation tank, separated by a 1 m high wall (the floc weir), which will have the advantage of eliminating the "lost triangle" from the plate settlers.
```python
TimeRun = VolStorage/FlowPlant # The floc hopper will be drained after every run.
FluxMaxSed = (10*u.NTU*1*u.L/u.s).to(u.kg/u.day) # Maximum mass of flocs being settled (conversion of NTU to mg/L assumes they are kaolin flocs, which is not likely the case.)
MassMaxHopper = FluxMaxSed*TimeRun
VolMinHopper = (MassMaxHopper/DensityCompactFloc).to(u.L)

ThickWeir = 10*u.cm
HeightWeir = 1*u.m
VolDesignHopper = ((0.2*u.m-ThickWeir)*HeightWeir*WidthSed).to(u.L)
VolMinHopper<VolDesignHopper # The proposed hopper volume seems to be more than adequate.
```
A hopper occupying the 0.2 m between the interior wall (floc weir) and the wall of the sedimentation tank seems to be adequate, even when considering a weir thickness of 10 cm occupying half of this volume.

### Effluent Manifold Design
```python
HeadLossLaunder = 4*u.cm # AguaClara constant
RoughnessPVC = 0.005*u.mm # Assumed, based on [Engineering Toolbox](https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html)
SpaceOrificeLaunderEstimated = 10*u.cm # AguaClara standard

DiamMinLaunder = (pc.manifold_id(FlowPlant,HeadLossLaunder,LengthSed,PiFlowManifold,nu,RoughnessPVC,1,np.floor((LengthSed/SpaceOrificeLaunderEstimated).to(u.dimensionless)))).to(u.inch)
LocPipeLaunder = (np.array(pipes.pipedb['ID_SCH40'])>DiamMinLaunder.magnitude).argmax() # Index for manifold in pipe_database.csv
DiamNomLaunder = pipes.pipedb.iloc[LocPipeLaunder,0]*u.inch
print('The nominal diameter of the exit launder is ',DiamNomLaunder,'.',sep='')

NumOrificeLaunderEstimated = (LengthSed/SpaceOrificeLaunderEstimated).to(u.dimensionless)

DiamOrificeLaunder = (ut.ceil_nearest(pc.diam_circle((FlowPlant/NumOrificeLaunderEstimated)/(con.VC_ORIFICE_RATIO*np.sqrt(2*u.g_0*HeadLossLaunder))),drills.DRILL_BITS_D_METRIC)).to(u.mm)

print('The diameter of the launder orifices is ',DiamOrificeLaunder,'.',sep='')

FlowOrificeLaunder = pc.flow_orifice_vert( DiamOrificeLaunder, HeadLossLaunder, con.VC_ORIFICE_RATIO )

SpaceOrificeLaunder = ((LengthSed - pipe.socket_depth(DiamNomLaunder) - pipe.cap_thickness(DiamNomLaunder) - DiamOrificeLaunder) / ((FlowPlant / FlowOrificeLaunder) - 1)).to(u.cm)
print('The space between orifices on the launder is ',np.round(SpaceOrificeLaunder,1),'.',sep='')

NumOrificeLaunder = math.floor((LengthSed - pipe.socket_depth(DiamNomLaunder) - pipe.cap_thickness(DiamNomLaunder) - DiamOrificeLaunder) / SpaceOrificeLaunder ) + 1
print('The number of orifices on the launder is ',NumOrificeLaunder,'.',sep='')

DiamNomScupper = 2*u.inch
print('The scupper diameter is ',DiamNomScupper,'.',sep='')
```
The nominal diameter of the exit launder is 2 inches with a length equal to the length of the flocculator. The orifices are 14 mm in diameter spaced apart at 10.8 cm for a total of 11 orifices. For a scupper (skimmer) located at the free surface of the tank, its diameter is assumed to be 2 inches.

### Height of Sedimentation Tank
```python
# HeightManifold = (0.5*pipes.pipedb.iloc[LocPipeManifold,3]*u.mm).to(u.cm)
# print('The height of the influent diffuser is ',np.round(HeightManifold,1),'.',sep='')

HeightSlope = (WidthSed/2)*np.tan(AngleSlope)
print('The height of the sloped bottoms is ',np.round(HeightSlope,2),'.',sep='')

HeightWeir2Plate = 10*u.cm # Assume 10 cm from top of floc weir to bottom of plate settlers
HeightPlate2Launder = 5*u.cm # Height from top of plate settlers to bottom of launder.
HeightFreeBoard = 5*u.cm
HeightSed = HeightWeir + HeightWeir2Plate + HeightPlate + HeightPlate2Launder + pipes.pipedb.iloc[LocPipeLaunder,1]*u.inch + HeightFreeBoard
print('The height of the sedimentation tank is ',np.round(HeightSed,2),'.',sep='')

HeightOutfall = 5*u.cm
HeightSumpTank = HeightWeir + HeightWeir2Plate + HeightPlate + HeightPlate2Launder - HeightOutfall
print('The height of the sump tank is ',np.round(HeightSumpTank,2),'.',sep='')
```
The sedimentation tank is calculated to have a total height of 1.59 m. Assuming that the storage tank will be 5 cm below the bottom of the sedimentation launder, the storage tank will be at a maximum 1.43 m tall.
