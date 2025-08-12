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
    """
    Plot air flow vs depth for different air pressures
    """
    depth_min, depth_max = depth_lims
    depths = np.linspace(*depth_lims)
    
    # Debug: Print to check inputs
    print(f"Debug - Depth range: {depth_min} to {depth_max}")
    print(f"Debug - Number of air pressures: {len(air_pressures)}")
    print(f"Debug - Air pressures: {air_pressures[:3] if len(air_pressures) > 3 else air_pressures}...")
    
    if units == 'SI':
        column_names = ['depth (m)', 'air flow (kg/s)', 'air flow (SCMM)', 'air pressure (kPa)']
        pressure_unit = 'kPa'
    elif units == 'English':
        column_names = ['depth (ft)', 'air flow (kg/s)', 'air flow (SCFM)', 'air pressure (psi)']
        pressure_unit = 'psi'
    
    # Initialize empty DataFrame with proper dtypes
    source1 = pd.DataFrame(columns=column_names)
    
    # Collect all dataframes first, then concatenate once
    dataframes_to_concat = []
    
    for i, air_pressure in enumerate(air_pressures):
        print(f"Debug - Processing pressure {i+1}/{len(air_pressures)}: {air_pressure}")
        
        # Your existing calculation code
        solved_geometries = [
            copy.deepcopy(
                downstream_solve(
                    system_geometry=system_geometry,
                    air_pressure=air_pressure,
                    water_pressure=convert.pressure_to_absolute(convert.H_m_to_Pa(depth)),
                    starting_mdot=0.1,
                    air_temp=air_temp
                )
            ) for depth in depths
        ]
        
        air_flows_mass = [parse_results(geom).total_mdot() for geom in solved_geometries]
        
        atm_pressure = convert.pressure_to_absolute(0)
        rho_air = air.rho_air(atm_pressure, T=convert.F_to_C(68))
        
        air_flows_volumetric = [convert.CMS_to_CMM(air.Q(mdot, rho_air)) for mdot in air_flows_mass]
        
        # Create dataframe for this pressure
        if units == 'SI':
            df_temp = pd.DataFrame({
                column_names[0]: depths,
                column_names[1]: air_flows_mass,
                column_names[2]: air_flows_volumetric,
                column_names[3]: [convert.pressure_to_gauge(air_pressure)/1000] * len(depths)  # Repeat for all depths
            })
        elif units == 'English':
            df_temp = pd.DataFrame({
                column_names[0]: list(map(convert.m_to_ft, depths)),
                column_names[1]: air_flows_mass,
                column_names[2]: list(map(convert.CMM_to_CFM, air_flows_volumetric)),
                column_names[3]: [convert.Pa_to_psi(convert.pressure_to_gauge(air_pressure))] * len(depths)  # Repeat for all depths
            })
        
        dataframes_to_concat.append(df_temp)
    
    # Concatenate all dataframes at once
    if dataframes_to_concat:
        source1 = pd.concat(dataframes_to_concat, ignore_index=True)
    
    # Debug: Check the resulting dataframe
    print(f"Debug - Final dataframe shape: {source1.shape}")
    print(f"Debug - Columns: {source1.columns.tolist()}")
    print(f"Debug - First few rows:\n{source1.head()}")
    print(f"Debug - Data types:\n{source1.dtypes}")
    print(f"Debug - Unique pressures: {source1[column_names[3]].unique()}")
    
    # Ensure numeric columns are actually numeric
    source1[column_names[0]] = pd.to_numeric(source1[column_names[0]], errors='coerce')
    source1[column_names[1]] = pd.to_numeric(source1[column_names[1]], errors='coerce')
    source1[column_names[2]] = pd.to_numeric(source1[column_names[2]], errors='coerce')
    source1[column_names[3]] = pd.to_numeric(source1[column_names[3]], errors='coerce')
    
    # Check for NaN values
    if source1.isnull().any().any():
        print("Warning: DataFrame contains NaN values!")
        print(source1.isnull().sum())
    
    # Create the base line chart
    chart = alt.Chart(source1).mark_line(point=True).encode(
        x=alt.X(column_names[0], title=column_names[0]),
        y=alt.Y(column_names[2], title=column_names[2]),
        color=alt.Color(
            column_names[3],
            title=f'Air Pressure ({pressure_unit})',
            scale=alt.Scale(scheme='viridis'),
            legend=alt.Legend(orient='right')
        ),
        tooltip=[
            alt.Tooltip(column_names[0], format='.2f'),
            alt.Tooltip(column_names[2], format='.2f'),
            alt.Tooltip(column_names[3], format='.1f')
        ]
    )
    
    # Create labels for each pressure line
    # Group by pressure and get the middle point for label placement
    label_data = source1.groupby(column_names[3]).agg({
        column_names[0]: 'median',
        column_names[2]: 'median'
    }).reset_index()
    
    # Add formatted label column
    label_data['label'] = label_data[column_names[3]].apply(lambda x: f'{x:.0f} {pressure_unit}')
    
    labels = alt.Chart(label_data).mark_text(
        align='left',
        dx=5,
        dy=-5,
        fontSize=12
    ).encode(
        x=column_names[0],
        y=column_names[2],
        text='label:N',
        color=alt.Color(
            column_names[3],
            scale=alt.Scale(scheme='viridis'),
            legend=None
        )
    )
    
    # Combine chart and labels
    final_chart = (chart + labels).properties(
        title='Air Flow vs Water Depth',
        height=600,
        width=800
    ).configure_axis(
        labelFontSize=12,
        titleFontSize=14
    ).configure_title(
        fontSize=16
    ).interactive()  # Add interactivity for zooming/panning
    
    return final_chart

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
                size=80,
                opacity=0.8
            ).encode(
            x = alt.X(orifice_label, scale=alt.Scale(domain=xdomain)), 
            y= alt.Y(flow_label, scale=alt.Scale(domain=ydomain))
            ).properties(title='Orifice Flows', height=800, width=800
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
                size=80,
                opacity=0.8
            ).encode(
            x = alt.X(orifice_label, scale=alt.Scale(domain=xdomain)),
            y= alt.Y(vel_label, scale=alt.Scale(domain=ydomain))
            ).properties(title='Surface Velocities', height=800, width=800
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