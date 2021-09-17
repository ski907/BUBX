from geometry import Pipe
from geometry import Orifice
from physics import air

import numpy as np

from conversions import convert



class parse_results:
    def __init__(self, system_geometry):
        self.system_geometry = system_geometry
        self.atm_pressure = convert.pressure_to_absolute(0)
        self.rho_air_standard = air.rho_air(self.atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
        
    def orifice_positions(self):
        offsets = [0]
        for i,p in enumerate(self.system_geometry.pipes):
            offsets.append(p.length + offsets[i]) 
        return offsets[1:]
        
    
    def total_mdot(self):
        return sum([orifice.mdot for orifice in self.system_geometry.orifices])
    
    def orifice_mdots(self):
        return [o.mdot for o in self.system_geometry.orifices]
    
    def orifice_flows_SCMM(self):
        return [convert.CMS_to_CMM(air.Q(mdot,self.rho_air_standard)) for mdot in self.orifice_mdots()]
    
    def orifice_flows_SCFM(self):
        return [convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot,self.rho_air_standard))) for mdot in self.orifice_mdots()]
    
    def orifice_flows_SCMM_per_unit_length(self):
        pipe_lengths = [p.length for p in self.system_geometry.pipes][1:]
        spacings = [(a + b) / 2 for a, b in zip(pipe_lengths[:], pipe_lengths[1:])]
        spacings.insert(0,pipe_lengths[0])
        spacings.append(pipe_lengths[-1])
        return [convert.CMS_to_CMM(air.Q(mdot,self.rho_air_standard))/spacings[i] for i, mdot in enumerate(self.orifice_mdots())]
    
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
        #modified from Haehnel 2016 to go to 1 when perfectly uniform
        #it's just 1- CU from Haehnel
        return 1 - (self.system_geometry.orifices[0].mdot - self.system_geometry.orifices[-1].mdot)/self.mean_orifice_mass_flow()
               
    def mean_orifice_mass_flow(self):
        return np.mean([o.mdot for o in self.system_geometry.orifices])
        
    def pipe_mach_numbers(self):
        return [p.mach_number for p in self.system_geometry.pipes]

    def pipe_reynolds_numbers(self):
        return [p.Re for p in self.system_geometry.pipes]
    
    def airflow_per_unit_length(self):
        mdot_per_unit_length =  self.total_mdot() / sum([p.length for p in self.system_geometry.pipes[1:]])
    
        return convert.CMS_to_CMM(air.Q(mdot_per_unit_length,self.rho_air_standard))
    
    def horizontal_surface_vel_haehnel2016(self):
        water_depth = convert.Pa_to_H_m(
                                        convert.pressure_to_gauge(
                                        self.system_geometry.orifices[0].downstream_pressure)
                                        )
        Qa = self.airflow_per_unit_length()
        
        return ( (351+25.47*water_depth) / 1000) * (Qa ** 0.4223)

    def horizontal_surface_vel_haehnel2016_at_orifices(self):
        water_depth = convert.Pa_to_H_m(
                                        convert.pressure_to_gauge(
                                        self.system_geometry.orifices[0].downstream_pressure)
                                        )
        Qa_at_orifices = self.orifice_flows_SCMM_per_unit_length()
        
        return [( (351+25.47*water_depth) / 1000) * (Qa ** 0.4223) for Qa in Qa_at_orifices]   
    
    
    def orifice_to_diffuser_area_ratio(self):
        total_orifice_area = sum([o.area for o in self.system_geometry.orifices])
        diffuser_pipe_area = self.system_geometry.pipes[-1].area
        ratio = total_orifice_area / diffuser_pipe_area
        
        return ratio
            