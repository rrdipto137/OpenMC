import openmc
import numpy as np
import math
from math import pi, sin, cos
from matplotlib import pyplot
import matplotlib.pyplot as plt






#material file

#fuel
uranium = openmc.Material()
uranium.add_nuclide('U235', 4.2)
uranium.add_nuclide('U238', 95.8)


o = openmc.Material()
o.add_element('O', 1.0)

uo2 = openmc.Material.mix_materials([o, uranium],[.6633, 0.3367])


uranium2 = openmc.Material()
uranium2.add_nuclide('U235', 3.7)
uranium2.add_nuclide('U238', 96.3)


LEU = openmc.Material.mix_materials([o, uranium2],[.6633, 0.3367])



#cladding
#AIM 1 Ti-15-15 steel

clad = openmc.Material()
clad.add_element('Cr', 15, 'wo')
clad.add_element('Ni', 15, 'wo')
clad.add_element('Mo', 1.5, 'wo')
clad.add_element('Mn', 1.5, 'wo')
clad.add_element('Si', 0.9, 'wo')
clad.add_element('Ti', 0.4, 'wo')
clad.add_element('C', 0.09, 'wo')
clad.add_element('B', 0.006, 'wo')
clad.add_element('Fe', 65.604,'wo')
clad.set_density('g/cm3', 7.95)
clad.temperature = 900



#gap
he = openmc.Material()
he.add_element('He', 1)
he.set_density('g/cm3', 0.00016)
he.temperature= 700
 
#absorber

gd_u1 = openmc.Material(name="UO2-5%Gd2O3, 3.3% enriched")
gd_u1.set_density('g/cm3', 10.4)
gd_u1.temperature = 700

# uranium
gd_u1.add_element('U', 12.75, enrichment=3.3)
# gadolinium
gd_u1.add_element('Gd', 1.0)
# oxygen
gd_u1.add_element('O', 27.0)


zirconium= openmc.Material()
zirconium.add_elements_from_formula('ZrO2','ao', 1.0)

yt = openmc.Material()
yt.add_elements_from_formula('Y2O3', 'ao', 1.0)

zirconia = openmc.Material.mix_materials([zirconium, yt],[0.95, 0.05], 'wo')
zirconia.temperature=700
zirconia.set_density('g/cm3', 6.08)

mat= openmc.Materials([uo2, LEU,clad, he,  gd_u1,  zirconia])
mat.export_to_xml()

#geometry

fuel_ir = openmc.ZCylinder( r=.1)
fuel_or1 = openmc.ZCylinder( r=.45) 
fuel_or2 = openmc.ZCylinder( r=.45) 

dummy = openmc.ZCylinder( r=.45) 


gd_or = openmc.ZCylinder(r=.45)

clad_ir1 = openmc.ZCylinder( r=.465) 
clad_or1 = openmc.ZCylinder( r=.525) 

clad_ir2 = openmc.ZCylinder( r=.465) 
clad_or2 = openmc.ZCylinder( r=.525) 



top = openmc.ZPlane(z0=+191.85)
bottom = openmc.ZPlane(z0=-191.85) 

up = openmc.ZPlane(z0=+195,boundary_type='vacuum')
down = openmc.ZPlane(z0=-195, boundary_type='vacuum')

gd_region= -gd_or  & -up & +down




fuel_region1 = -fuel_or1 & -up & +down
gap_region1  = +fuel_or1 & -clad_ir1  & -top & +bottom
clad_region = +clad_ir1 & -clad_or1  & -top & +bottom
coolant_region1 = +clad_or1 & -up & +down

coolant_gd = +gd_or & -up & +down
dummy_region= -dummy & -up & +down

fuel_region2 = -fuel_or2 & -up & +down
gap_region2  = +fuel_or2 & -clad_ir2  & -top & +bottom
clad_region2 = +clad_ir2 & -clad_or2  & -top & +bottom
coolant_region2 = +clad_or2 & -up & +down


coolant_dummy = + dummy & -up & +down



#MOX1
fuel_cell1 = openmc.Cell(name="fuel_cell1", fill=uo2, region=fuel_region1)
gd_cell1 = openmc.Cell(name="gd_cell1", fill=gd_u1, region=gd_region)
gap_cell1 = openmc.Cell(name="gap_cell1",fill=he, region=gap_region1)
clad_cell1 = openmc.Cell(name="clad_cell1", fill=clad, region=clad_region)
coolant_cell1 = openmc.Cell(name="coolant_cell1", fill=he, region=coolant_region1)

coolant_gd = openmc.Cell( name="lk_cell1",fill=he, region=coolant_gd)

dummy_cell1 = openmc.Cell(name="sd_cell1", fill=he, region=dummy_region)
coolant_dummy_cell = openmc.Cell(name="gsdcell1", fill=he, region=coolant_dummy)

#MOX2
fuel_cell2 = openmc.Cell(name="gddsell1", fill=LEU, region=fuel_region2)
gap_cell2 = openmc.Cell(name="gdadsll1",fill=he, region=gap_region2)
clad_cell2 = openmc.Cell(name="gd1", fill=clad, region=clad_region2)
coolant_cell2 = openmc.Cell(name="gdadsll1", fill=he, region=coolant_region2)




fuel_u1 = openmc.Universe(cells=[fuel_cell1, gap_cell1, clad_cell1, coolant_cell1])
fuel_u2 = openmc.Universe(cells=[fuel_cell2, gap_cell2, clad_cell2, coolant_cell2])
dummy_l = openmc.Universe(cells=[gd_cell1,coolant_gd])
dummyp = openmc.Universe(cells=[dummy_cell1,coolant_dummy_cell])
'''
fuel_u1.plot(width=(1,1), pixels = (1000,1000), color_by='material',colors= {uo2 : 'orange', he : 'blue',LEU : 'pink',gd_u1 : 'black'})
fuel_u2.plot(width=(1,1), pixels = (1000,1000), color_by='material',colors= {uo2 : 'orange', he : 'blue',LEU : 'pink',gd_u1 : 'black'})
dummy_l.plot(width=(1,1), pixels = (1000,1000), color_by='material',colors= {uo2 : 'orange', he : 'blue',LEU : 'pink',gd_u1 : 'black'})

dummyp.plot(width=(1,1), pixels = (1000,1000), color_by='material',colors= {uo2 : 'orange', he : 'blue',LEU : 'pink',gd_u1 : 'black'})
'''
#inner fuel assembly
all_lead_out = openmc.Cell(fill=he)
all_lead_out_u = openmc.Universe(cells=[all_lead_out])
in_lat=openmc.HexLattice(name='assembly')
in_lat.center = (0., 0.)
in_lat.pitch = (1.275,)
in_lat.outer=all_lead_out_u

in_ring11 = [fuel_u1]*66
in_ring10 = [fuel_u2]*60
in_ring9 = [fuel_u2,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,fuel_u2,fuel_u1,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u2,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u2,
            fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,fuel_u2,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,fuel_u1,fuel_u2,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1
            ]

in_ring8 = [fuel_u1]*48

in_ring7 = [fuel_u1]*42

in_ring6 = [fuel_u1, fuel_u1, fuel_u1,      
  dummyp, fuel_u1, fuel_u1,      
    fuel_u1, fuel_u1, fuel_u1,      
  dummyp, fuel_u1, fuel_u1, 
  fuel_u1, fuel_u1, fuel_u1,
   dummyp, fuel_u1, fuel_u1,
    fuel_u1, fuel_u1, fuel_u1,
  dummyp,fuel_u1, fuel_u1,
    fuel_u1, fuel_u1, fuel_u1,
  dummyp,fuel_u1, fuel_u1,
    fuel_u1, fuel_u1, fuel_u1,
   dummyp,  fuel_u1, fuel_u1]

in_ring5 = [dummyp,fuel_u1,fuel_u1,fuel_u1,fuel_u1,dummyp,
            fuel_u1,fuel_u1,fuel_u1,fuel_u1,dummyp,fuel_u1,
            fuel_u1,fuel_u1,fuel_u1,dummyp,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,dummyp,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,dummyp,fuel_u1,fuel_u1,fuel_u1,fuel_u1
    
    
    
    
    ]

in_ring4 = [ fuel_u1,fuel_u1,fuel_u1,dummy_l,fuel_u1,fuel_u1,
            fuel_u1,dummy_l,fuel_u1,fuel_u1,fuel_u1,dummy_l,
            fuel_u1,fuel_u1,fuel_u1,dummy_l,fuel_u1,fuel_u1,
            fuel_u1,dummy_l,fuel_u1,fuel_u1,fuel_u1,dummy_l]

in_ring3 = [fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,fuel_u1,
            fuel_u1,fuel_u1,fuel_u1,fuel_u1 ]
in_ring2 = [fuel_u1]*12
in_ring1 = [fuel_u1]*6
in_ring0 = [dummyp]

in_lat.universes = [in_ring11,in_ring10, in_ring9, in_ring8, in_ring7,in_ring6,in_ring5,in_ring4, in_ring3, in_ring2, in_ring1, in_ring0]
in_lat.orientation='y'


outer_in_surface = openmc.model.HexagonalPrism(edge_length=15, orientation='y')

main_in_assembly = openmc.Cell(fill=in_lat, region=-outer_in_surface & -top & +bottom)



in_plenum_top = openmc.Cell(fill= he, region= -outer_in_surface & -up & +top )
in_plenum_bottom = openmc.Cell(fill= he, region= -outer_in_surface & +down & -bottom )

out_in_assembly  = openmc.Cell(fill=he, region=+outer_in_surface & -up & +down)

main_in_u = openmc.Universe( cells=[main_in_assembly, out_in_assembly, in_plenum_bottom, in_plenum_top])

main_in_u.plot(width=(30,30), pixels = (1000,1000), color_by='material',colors= {uo2 : 'orange', he : 'blue',LEU : 'pink',gd_u1 : 'black'})
'''

outer_empty_surface = openmc.model.hexagonal_prism(edge_length=9, orientation='y')
empty_cell= openmc.Cell(fill= he, region= outer_empty_surface & -up & +down)
out_empty_cell  = openmc.Cell(fill=he,region= ~outer_empty_surface & -up & +down)
empty_u= openmc.Universe(cells=[empty_cell,out_empty_cell])


core_lat = openmc.HexLattice( name='core')
core_lat.center = (0., 0.)
core_lat.pitch = (28.6,)
core_lat.outer = all_lead_out_u
core_lat.orientation ='x'

core_ring5 = [empty_u]*30
core_ring4 = [empty_u]*24
core_ring3 = [empty_u]*10 + [main_in_u]*8 
core_ring2 = [main_in_u]*2 + [empty_u]*3 + [main_in_u]*7 
core_ring1 = [main_in_u]*6
core_ring0 = [main_in_u]*1
core_lat.universes = [core_ring5, core_ring4, core_ring3, core_ring2, core_ring1, core_ring0 ]

outer_core_surface = openmc.Cylinder(x0=0, y0=19, r=80, boundary_type='vacuum')

core = openmc.Cell( fill=core_lat, region=-outer_core_surface & -up & +down)

main_u = openmc.Universe( cells=[core]) 
main_u.plot(width= (160,160), pixels = (800,800),color_by='cell')
main_u.plot(width=(400,400), pixels=(500,500), color_by='cell',basis='xz')
main_u.plot(width=(400,400), pixels=(500,500), color_by='cell',basis='xy')


'''
'''
geom= openmc.Geometry(main_u)
geom.export_to_xml()



settings = openmc.Settings()

settings.temperature = {'method': 'interpolation'}
settings.batches = 150
settings.inactive = 5
settings.particles =900
settings.export_to_xml()




model = openmc.Model(materials=mat,geometry=geom, settings=settings)


d = 700 
# Instantiate a tally Mesh
mesh = openmc.RegularMesh()
mesh.dimension = [d, d]
mesh.lower_left = [-80, -80]
mesh.upper_right = [+80, +80]

# Instantiate tally Filter
mesh_filter = openmc.MeshFilter(mesh)

# Creating Conversion Factor of Flux
# Creating a tally of fission-q-recoverable
H = openmc.Tally()
H.scores = ['heating-local']
H.filter = [mesh_filter]

# Instantiate the Tally
tally = openmc.Tally(name='mesh tally')
tally.filters = [mesh_filter]
tally.scores = ['flux']
#tally.scores = ['flux','fission']     # takes fission and flux
tallies = openmc.Tallies([H, tally])

model.tallies = tallies

# Run OpenMC
statepoint_filename = model.run()

ce_spfile = './statepoint.150.h5'
#os.rename(statepoint_filename, ce_spfile)
# Move the Summary file
ce_sumfile = './summary.h5'
#os.rename('summary.h5', ce_sumfile)

# Load the statepoint file
sp = openmc.StatePoint(ce_spfile, autolink=False)

H_value = sp.tallies[H.id].get_pandas_dataframe()['mean'][0]

# ================================================================
# Calculating the Volume of Mesh
mesh_height = 60
mesh_area = (160 / d)**2

volume = mesh_height * mesh_area

# Calculating the Factor to Normalize the Flux
power = 300e6
H_Jsrc = 1.602e-19 * H_value

f = power / H_Jsrc
flux_normalization_factor = f/volume

# ================================================================



# Load the summary file in its new location
su = openmc.Summary(ce_sumfile)
sp.link_with_summary(su)

# Get the OpenMC fission rate mesh tally data
ce_mesh_tally = sp.get_tally(name='mesh tally')
ce_flux = ce_mesh_tally.get_values(scores=['flux'])

# Reshape array to 2D for plotting
ce_flux.shape = mesh.dimension
print(ce_flux)
# Normalize to the average pin power
ce_flux /= np.mean(ce_flux[ce_flux > 0.])
print(ce_flux)
# Force zeros to be NaNs so their values are not included when matplotlib calculates
# the color scale
ce_flux[ce_flux == 0.] = np.nan


plt.figure(dpi=1200)
plt.imshow(ce_flux*flux_normalization_factor, interpolation='none', origin='lower', cmap = 'twilight_shifted')
plt.colorbar()
plt.title('Flux Distribution')
plt.show()
'''