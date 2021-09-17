from solver import downstream_solve

import pandas as pd
import numpy as np
import altair as alt
from physics import air

from conversions import convert

from geometry.segments import generate_equal_spaced_uniform
from geometry.parse_results import parse_results

import copy

def plot_flow_vs_depth(system_geometry,
                       depth_lims: tuple,
                       air_pressures: list,
                       air_temp,
                       units='SI'):

    depth_min, depth_max = depth_lims
    depths = np.linspace(*depth_lims)
    
    if units == 'SI':
        column_names = ['depth (m)', 'air flow (kg/s)', 'air flow (SCMM)', 'air pressure (kPa)']
        pressure_unit = 'kPa'
    elif units == 'English':
        column_names = ['depth (ft)', 'air flow (kg/s)', 'air flow (SCFM)', 'air pressure (psi)']
        pressure_unit = 'psi'
        
    source1 = pd.DataFrame(columns = column_names, dtype=object)

    for air_pressure in air_pressures:    
#        air_flows_mass = [get_just_airflow(system_geometry = system_geometry, 
#                         air_pressure =  air_pressure, 
#                         water_pressure = convert.pressure_to_absolute(convert.H_m_to_Pa(depth)),
#                         starting_mdot= 0.1, 
#                         air_temp = air_temp) for depth in depths]
        
        solved_geometries = [copy.deepcopy(downstream_solve(system_geometry = system_geometry, 
                     air_pressure =  air_pressure, 
                     water_pressure = convert.pressure_to_absolute(convert.H_m_to_Pa(depth)),
                     starting_mdot= 0.1, 
                     air_temp = air_temp)) for depth in depths]
    
        air_flows_mass = [parse_results(geom).total_mdot() for geom in solved_geometries]
        
        atm_pressure = convert.pressure_to_absolute(0)
        rho_air = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
        
        air_flows_volumetric = [convert.CMS_to_CMM(air.Q(mdot,rho_air)) for mdot in air_flows_mass]
        
        if units == 'SI':
            source1  = source1.append(pd.DataFrame({
                    column_names[0]: depths,
                    column_names[1]: air_flows_mass,
                    column_names[2]: air_flows_volumetric,
                    column_names[3]: convert.pressure_to_gauge(air_pressure)/1000
                    }))    
    
        if units == 'English':
            source1  = source1.append(pd.DataFrame({
                    column_names[0]: list(map(convert.m_to_ft,depths)),
                    column_names[1]: air_flows_mass,
                    column_names[2]: list(map(convert.CMM_to_CFM,air_flows_volumetric)),
                    column_names[3]: convert.Pa_to_psi(convert.pressure_to_gauge(air_pressure))
                    }))             
                
    chart = alt.Chart(source1).mark_line().encode(
            x=column_names[0],
            y=column_names[2],
            color=column_names[3]
            )
    
    labels = alt.Chart(source1).mark_text(align='left', dx=0, dy=-15).encode(
            alt.X(column_names[0], aggregate='mean', axis=alt.Axis(title=column_names[0])),
            alt.Y(column_names[2], aggregate='mean', axis=alt.Axis(title=column_names[2])),
            text = 'label:N',
            color = alt.Color(column_names[3], legend=None, scale=alt.Scale(domain=air_pressures,type='ordinal'))
            ).transform_calculate(label=f'format(datum["{column_names[3]}"],".0f") + " {pressure_unit}"')

    
    return alt.layer(chart, labels).properties(title='Range', height=800, width=1000).configure_axis(
                    labelFontSize=15,
                    titleFontSize=20
                    )

def plot_orifice_flows(solved_geom,units):
    
    if units == 'SI':
        flow_label = 'Orifice Airflow (SCMM)'
        orifice_label = 'orifice offset (m)'
        flows = parse_results(solved_geom).orifice_flows_SCMM()
        offsets = parse_results(solved_geom).orifice_positions()
    if units == 'English':
        flow_label = 'Orifice Airflow (SCFM)'
        orifice_label = 'orifice offset (ft)'
        flows = parse_results(solved_geom).orifice_flows_SCFM()
        offsets = [convert.m_to_ft(offset) for offset in parse_results(solved_geom).orifice_positions()]
    
    
    
    data = pd.DataFrame({flow_label:flows, orifice_label:offsets})
    
    xdomain = [min(data[orifice_label])*.95,max(data[orifice_label])*1.05]
    ydomain = [min(data[flow_label])*.95,max(data[flow_label])*1.05]
    chart = alt.Chart(data).mark_circle(
                color='red',
                opacity=0.8
            ).encode(
            x = alt.X(orifice_label, scale=alt.Scale(domain=xdomain)), 
            y= alt.Y(flow_label, scale=alt.Scale(domain=ydomain))
            ).properties(title='Orifice Flows', height=800, width=1000
            ).configure_axis(
                    labelFontSize=20,
                    titleFontSize=20
                    )
    
    return chart

def plot_horizontal_velocities(solved_geom,units):
    
    if units == 'SI':
        vel_label = 'Horizontal Surface Velocity (m/s)'
        orifice_label = 'orifice offset (m)'
        h_velocities = parse_results(solved_geom).horizontal_surface_vel_haehnel2016_at_orifices()
        offsets = parse_results(solved_geom).orifice_positions()
    if units == 'English':
        vel_label = 'Horizontal Surface Velocity (ft/s)'
        orifice_label = 'orifice offset (ft)'
        h_velocities = [convert.m_to_ft(hvel) for hvel in parse_results(solved_geom).horizontal_surface_vel_haehnel2016_at_orifices()]
        offsets = [convert.m_to_ft(offset) for offset in parse_results(solved_geom).orifice_positions()]
    
    
    
    data = pd.DataFrame({vel_label:h_velocities, orifice_label:offsets})
    
    xdomain = [min(data[orifice_label])*.95,max(data[orifice_label])*1.05]
    ydomain = [min(data[vel_label])*.95,max(data[vel_label])*1.05]
    chart = alt.Chart(data).mark_circle(
                color='red',
                opacity=0.8
            ).encode(
            x = alt.X(orifice_label, scale=alt.Scale(domain=xdomain)),
            y= alt.Y(vel_label, scale=alt.Scale(domain=ydomain))
            ).properties(title='Surface Velocities', height=800, width=1000
            ).configure_axis(
                    labelFontSize=20,
                    titleFontSize=20
                    )
    
    return chart
    
    
def main():
    
    manifold = generate_equal_spaced_uniform(segment_length=2800, 
                              pipe_roughness=0.000025,
                              pipe_diameter= 0.0762009266,
                              number_of_orifices=2,
                              orifice_diameter=0.015875193
                              )
    
    air_pressures = [30, 40, 50, 60, 70, 80, 90]
    
    air_pressures_Pa = [convert.pressure_to_absolute(convert.psi_to_Pa(pressure)) for pressure in air_pressures]
    plot_flow_vs_depth(manifold, (0,20), air_pressures_Pa, air_temp=0, units='English')

if __name__ == '__main__':
    main()