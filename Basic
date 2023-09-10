# Pincell-modeling-

import openmc
import matplotlib.pyplot as plt

# Materials definitions




# Fuel
fuel=openmc.Material(name='fuel' ,temperature=900)
fuel.add_element('U', 1.0, enrichment=4.2)
fuel.add_element('O', 2.0)
fuel.set_density('g/cc', 10.55)

# Moderator
water = openmc.Material(name='water',temperature=573)
water.add_nuclide('H1', 2.0)
water .add_nuclide('O16', 1.0)
water .set_density('g/cm3', 1.005)
water.add_s_alpha_beta('c_H_in_H2O')

# Cladding
zircalloy4 = openmc.Material(name='zircalloy4',temperature=613)
zircalloy4.add_element('Sn', 0.014, 'wo')
zircalloy4.add_element('Fe', 0.00165, 'wo')
zircalloy4.add_element('Cr', 0.001, 'wo')
zircalloy4.add_element('Zr', 0.98335, 'wo')
zircalloy4.set_density('g/cm3', 6.55)


# Gap_gas
helium=openmc.Material(name='helium',temperature=613)
helium.add_element('He', 1.0)
helium.set_density('g/cm3', 0.001598)

#hole_air

air=openmc.Material(name='air',temperature=900)
air.add_element('Ar', 0.01288)
air.add_element('N', 0.7552)
air.add_element('O', 0.2314)
air.add_element('C', 0.00051)
air.set_density('g/cm3', 0.0003922)

# Material File Generation
materials = openmc.Materials([fuel, water,zircalloy4, helium,air])
materials.export_to_xml()

# Geometry Definitions

bottom = openmc.ZPlane(z0=-1.0, boundary_type = 'reflective')
top = openmc.ZPlane(z0=1.0, boundary_type = 'reflective')
pitch = 1.275 #cm
pin_cell_box = openmc.rectangular_prism(width=pitch, height=pitch, boundary_type='reflective')

# Geometry definitions for the fuel rod
hole     = openmc.ZCylinder(r=0.15/2,name='hole')
fuel_or  = openmc.ZCylinder(r=0.757/2, name='Fuel OR')
fclad_ir = openmc.ZCylinder(r=(0.773/2), name='Clad IR')
fclad_or = openmc.ZCylinder(r=0.91/2, name='Clad OR')

hole_region=-hole
fuel_region = -fuel_or & +hole
gap_region  = +fuel_or & -fclad_ir
fclad_region  = +fclad_ir & -fclad_or
fwater_region = pin_cell_box & +fclad_or

hole_cell=openmc.Cell(name='hole')
hole_cell.fill=air
hole_cell.region=hole_region

fuel_cell = openmc.Cell(name='fuel')
fuel_cell.fill = fuel
fuel_cell.region = fuel_region 

gap_cell = openmc.Cell(name='gap')
gap_cell.fill = helium
gap_cell.region = gap_region

clad_cell = openmc.Cell(name='clad')
clad_cell.fill = zircalloy4 
clad_cell.region = fclad_region

fwater_cell = openmc.Cell(name='mod')
fwater_cell.fill = water
fwater_cell.region = fwater_region

fuel_pin_universe = openmc.Universe(cells=(fuel_cell,gap_cell,hole_cell, clad_cell, fwater_cell))

# Geometry File Generation
geometry = openmc.Geometry(fuel_pin_universe)
geometry.export_to_xml()

# Plotting
#plot = openmc.Plot.from_geometry(geometry)
#plot.color_by = 'cell'
#plot.to_ipython_image()

plot = openmc.Plot.from_geometry(geometry)
plot.color_by = 'material'
plot.colors = colors = {
    water: 'blue',
    fuel: 'yellow',
    air: 'red',
    zircalloy4 :'pink',
    helium:'black'
    
}
plot.to_ipython_image()

plot.pixels = (1024, 1024)
plot.to_ipython_image()


# OpenMC simulation parameters

point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)
settings = openmc.Settings()
settings.temperature={'method':'interpolation'}
settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 10000
settings.export_to_xml()

openmc.run()
