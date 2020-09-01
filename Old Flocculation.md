# Flocculator Design for Bhalunka

## Imports
```python
from aguaclara.play import *
import doctest
import aguaclara.core.pipes as pipes
import aguaclara.core.head_loss as hl
```

## General Parameters
```python
FlowPlant = 1*(u.L/u.s)
Temperature = 15*u.degC
nu = pc.viscosity_kinematic_water(Temperature)
```
## Original Design
This is a calculation using dimensions from the original [PF-300 fabrication manual](https://docs.google.com/document/d/1Lm7n5tZeKK7vLFa9XLA5ub1SRND5AQc6E5wPqYm3ZVw/edit). Reverse engineering the initial design seems to imply that it was designed for a Gθ of 37,000, with a head loss of about 2 m and a residence time of around 80 s.
### Nominal Velocity
```python
DiamFloc = (pipes.Pipe(nd=3*u.inch,sdr=26).id_sch40).to(u.cm) # Inner diameter of flocculator
AreaFloc = pc.area_circle(DiamFloc) # Area of pipe
VelocityFloc = (FlowPlant/AreaFloc).to(u.cm/u.s) # Nominal velocity through flocculator (without baffles)
print('The velocity through the flocculator is ',VelocityFloc,'.',sep='')
ReFloc = pc.re_pipe(FlowPlant,DiamFloc,nu)
print('The Reynolds number of the flocculator is ',ReFloc,'.',sep='')
VolumeFloc = AreaFloc*(227*u.cm + 8*170*u.cm + 16*6*u.cm) # Volume of the flocculator based on 1 inlet pipe, 8 channel pipes, and 8 pipe sections connecting the tees (rough approximation accounting for pipe elbows and pipe sections connecting tees as the same).
ThetaFloc = (VolumeFloc/FlowPlant).to(u.s) # Hydraulic residence time of the flocculator.
print('The hydraulic residence time of the flocculator is ',ThetaFloc,'.',sep='')
```

### Geometric Baffle Calculations
Because the baffles are circles with flat sections removed (cut by a secant), the dimensions of the baffles are calculated using [circular segment geometry](https://en.wikipedia.org/wiki/Circular_segment), which sets up a trigonometric calculation of the angle defining the resulting chord.
```python
DiamBaffleFloc = DiamFloc - 2*u.mm # Diameter of the baffles
AreaBaffleFloc = pc.area_circle(DiamPlateFloc) # Area of baffle (without flat section cut out)
RadiusPlateFloc = DiamPlateFloc/2 # Hypotenuse of the triangle defined by half the chord, the radius of the baffle, and the perpendicular bisector of the chord.
RadiusFlatPlateFloc = 5.5*u.cm - 4.05*u.cm # Adjacent of the previously mentioned triangle
HalfAngleFlatPlateFloc = np.cos(RadiusFlatPlateFloc/RadiusPlateFloc) # Half Angle calculated from the previously given dimensions.
AngleFlatPlateFloc = (HalfAngleFlatPlateFloc*2).to(u.deg) # Whole angle described by the radii of the circle and the chord.
print('The angle of the chord in the circular baffles is ',AngleFlatPlateFloc,'.',sep='')
```
### Estimate of Head Loss$
The calculation of a K value for the circular segment opening is approximate because it uses the [sudden expansion of a pipe](https://en.wikipedia.org/wiki/Borda%E2%80%93Carnot_equation) as a model.
```python
AreaConstrictFloc = RadiusPlateFloc**2/2*(AngleFlatPlateFloc-np.sin(AngleFlatPlateFloc)).to(u.rad) # Constricted area for pipe expansion calculation
DiamNomConstrictFloc = pc.diam_circle(AreaConstrictFloc) # Converting constricted area into a corresponding diameter, assuming that the area is circular for the calculation.
KPlateFloc = hl.k_value_orifice(pipe_id=DiamFloc,orifice_id=DiamNomConstrictFloc,orifice_l=0.25*u.inch,q=FlowPlant,nu=nu) # Calculation of K for a single baffle.
print('The K value of a single baffle is ',KPlateFloc,'.',sep='')
NumExpansionFloc = 27 # From design manual
HeadLossFloc = (NumExpansionFloc*KPlateFloc*VelocityFloc**2/2/u.g_0).to(u.m) # Calculation of total design head loss for original flocculator.
print('The design head loss for the flocculator is ',HeadLossFloc,'.',sep='')
```
### Calculation of $G$
The calculations show a $G$ of 456/s and a $G\theta$ of 36,600.
```python
EDRFloc = (u.g_0*HeadLossFloc/ThetaFloc).to(u.mW/u.kg) # Energy dissipation rate
print('The energy dissipation rate in the flocculator is ',EDRFloc,'.',sep='')
GFloc = (np.sqrt(EDRFloc/nu)).to(1/u.s) # Hydraulic velocity gradient
print('The hydraulic velocity gradient in the flocculator is ',GFloc,'.',sep='')
GThetaFloc = (GFloc*ThetaFloc).to(u.dimensionless) # Camp number
print('The Camp number for the flocculator is ', GThetaFloc,'.',sep='')
```
