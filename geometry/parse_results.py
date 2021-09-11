from geometry import Pipe
from geometry import Orifice
from physics import air

from conversions import convert

class parse_results:
    def __init__(self, system_geometry):
        self.system_geometry = system_geometry
        
    def total_mdot(self):
        return sum([orifice.mdot for orifice in self.system_geometry.orifices])
    
    def orifice_mdots(self):
        return [o.mdot for o in self.system_geometry.orifices]
    
    def pipe_velocities(self):
        return [p.velocity(p.Q) for p in self.system_geometry.pipes]
    
    def pipe_pressure_drops(self):
        return [p.pressure_drop for p in self.system_geometry.pipes]
    
    def pipe_total_pressure_drop(self):
        return sum([p.pressure_drop for p in self.system_geometry.pipes])
    
    def ratio_pressure_drop_to_input(self):
        
        upstream_feature = self.system_geometry.features[0]
        
        if isinstance(upstream_feature, Orifice):
            input_air_pressure = upstream_feature.upstream_pressure
        if isinstance(upstream_feature, Pipe):
            input_air_pressure = upstream_feature.p1
        
        return self.pipe_total_pressure_drop()/convert.pressure_to_gauge(input_air_pressure)
    
    def coefficient_of_uniformity(self):
        return self.system_geometry.orifices[-1].mdot\
               /self.system_geometry.orifices[0].mdot
        
    def pipe_mach_numbers(self):
        return [p.mach_number for p in self.system_geometry.pipes]

    def pipe_reynolds_numbers(self):
        return [p.Re for p in self.system_geometry.pipes]
    
    def airflow_per_unit_length(self):
        mdot_per_unit_length =  self.total_mdot() / sum([p.length for p in self.system_geometry.pipes])
    
        atm_pressure = convert.pressure_to_absolute(0)
        rho_air = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
        
        return convert.CMS_to_CMM(air.Q(mdot_per_unit_length,rho_air))