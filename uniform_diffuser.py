import streamlit as st

from physics import air
from conversions import convert

from geometry.segments import generate_equal_spaced_uniform
from geometry.segments import add_supply_line

from solver import downstream_solve
from geometry.parse_results import parse_results

def uniform_diffuser():
    

    
    col1, col2, col3 = st.beta_columns(3)
    
    col1.write('System Geometry')
    supply_choice = col1.checkbox(label='Include Supply Line?')
    
    with col1.form(key='system_geometry_form'):
    
        st.write('Diffuser Geometry Input:')
        #pipe_roughness = st.number_input(label='Pipe Roughness Height (m)', value=0.000025, format='%1.10f')
        
        pipe_materials = {'stainless steel: 0.015mm ':0.015,
                          'galvinized steel: 0.15mm':0.15,
                          'PVC: 0.0015mm':0.0015,
                          'bub300 Stainless Steel: 0.025':0.025,
                          'Wire Reinforced Rubber Hose: 0.3mm':0.3}
        pipe_type = st.selectbox(label='Pipe Inner Material (Roughness in mm)',options=list(pipe_materials.keys()))
        pipe_roughness = pipe_materials[pipe_type]/1000 #mm to meters
        
        pipe_diameter = st.number_input(label='Pipe Diameter (in)', value=3., format='%1.2f')
        pipe_diameter = convert.in_to_m(pipe_diameter)
        
        segment_length = st.number_input(label='Segment Length (ft)', value=100)
        segment_length = convert.ft_to_m(segment_length)
        
        number_of_orifices = st.number_input(label='Number of Orifices', value=10, format='%d')
        
        orifice_diameter = st.number_input(label='Orifice Diameter (in)', value=5/8, format='%1.4f')
        orifice_diameter = convert.in_to_m(orifice_diameter)
        
        if supply_choice:
            st.write('Supply Line Geometry Input:')
            supply_pipe_type = st.selectbox(label='Supply Pipe Inner Material',options=list(pipe_materials.keys()))
            supply_pipe_roughness = pipe_materials[supply_pipe_type]/1000 #mm to meters
        
            supply_pipe_diameter = st.number_input(label='Supply Pipe Diameter (in)', value=3., format='%1.2f')
            supply_pipe_diameter = convert.in_to_m(supply_pipe_diameter)
        
            supply_pipe_length = st.number_input(label='Supply Pipe Length (ft)', value=100)
            supply_pipe_length = convert.ft_to_m(supply_pipe_length)
        
        sg_submit_button = st.form_submit_button(label='Update/Run')
    
    # st.form_submit_button returns True upon form submit
    if sg_submit_button:
        col1.write(f'Geometry Updated')
        
    manifold = generate_equal_spaced_uniform(segment_length, 
                                  pipe_roughness,
                                  pipe_diameter,
                                  number_of_orifices,
                                  orifice_diameter
                                  )
    
    if supply_choice:
        manifold =  add_supply_line(manifold,
                                    supply_pipe_length,
                                    supply_pipe_roughness,
                                    supply_pipe_diameter)
    
    col2.write('Boundary Conditions')
    
    with col2.form(key='Boundary Conditions'):
        air_pressure = st.number_input(label='Air Pressure (psi)', value=90, format='%d')
        air_pressure = convert.psi_to_Pa(air_pressure)
        air_pressure = convert.pressure_to_absolute(air_pressure)
        
        water_depth = st.number_input(label='Water Depth (ft)', value=30.0, format='%1.1f')
        water_depth = convert.ft_to_m(water_depth)
        water_pressure = convert.pressure_to_absolute(convert.H_m_to_Pa(water_depth))
        
        air_temp = st.number_input(label='Air Temp (C)', value=0, format='%d')
        
        starting_mdot = st.number_input(label='Starting Airflow (kg/s)', value=0.01)
        
        bc_submit_button = st.form_submit_button(label='Update/Run')
    
    # st.form_submit_button returns True upon form submit
    if bc_submit_button:
        col2.write(f'Boundary Conditions Updated')
    
    
    solved_geom = downstream_solve(manifold,
                           air_pressure, 
                           water_pressure, 
                           starting_mdot, 
                           air_temp,
                           verbose=False)
    
    #mdot = sum([orifice.mdot for orifice in solved_geom.orifices])
    
    col3.write('Results:')
    mdot = parse_results(solved_geom).total_mdot()
    col3.write('Total Mass flow rate: {:.2f} kg/s'.format(mdot))
    
    atm_pressure = convert.pressure_to_absolute(0)
    rho_air = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
    
    col3.write('Total Flow Rate (SI): {:.2f} SCMS @1atm and 20C'.format(air.Q(mdot,rho_air)))
    
    col3.write('Total Flow Rate (IMP): {:,.2f} SCFM @1atm and 68F'.format(convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot,rho_air)))))
    
    rho_air = air.rho_air(atm_pressure, T=0) #0C is Standard Temperature in SCFM in BUB300...
    col3.write('Total Flow Rate (IMP): {:,.2f} SCFM @1atm and 32F'.format(convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot,rho_air)))))

    col3.write('Airflow Per Unit Length: {:,.2f} SCMM/m'.format(parse_results(solved_geom).airflow_per_unit_length()))
    show_detailed_output = col3.checkbox(label='Show Detailed Output?')
    
    if show_detailed_output:
        col3.write('Total Pipe Pressure Drop (psi)')
        total_p_drop = convert.Pa_to_psi(parse_results(solved_geom).pipe_total_pressure_drop())
        col3.write(total_p_drop)
        
        col3.write('Ratio of Total Pressure Drop to Input Pressure')
        col3.write(parse_results(solved_geom).ratio_pressure_drop_to_input())
        
        col3.write('Coefficient of Uniformity')
        col3.write(parse_results(solved_geom).coefficient_of_uniformity())
        
        col3.write('Mass Flow Per Orifice (kg/s)')
        col3.write(parse_results(solved_geom).orifice_mdots())
        
        col3.write('Pipe Segment Velocities (m/s)')
        col3.write(parse_results(solved_geom).pipe_velocities())
        
        col3.write('Pipe Segment Pressure Drops (psi)')
        col3.write(list(map(convert.Pa_to_psi,parse_results(solved_geom).pipe_pressure_drops())))
        
        col3.write('Pipe Segment Mach Numbers')
        col3.write(parse_results(solved_geom).pipe_mach_numbers())
        
        col3.write('Pipe Segment Reynolds Number')
        col3.write(parse_results(solved_geom).pipe_reynolds_numbers())
    
    return manifold