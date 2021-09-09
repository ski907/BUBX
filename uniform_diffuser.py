import streamlit as st

from physics import air
from conversions import convert

from geometry.segments import generate_equal_spaced_uniform

from solver import downstream_solve

def uniform_diffuser():

    col1, col2 = st.beta_columns(2)
    
    col1.write('System Geometry')
    
    with col1.form(key='system_geometry_form'):
    
        
        #pipe_roughness = st.number_input(label='Pipe Roughness Height (m)', value=0.000025, format='%1.10f')
        
        pipe_materials = {'stainless steel':0.015,
                          'galvinized steel':0.15,
                          'PVC':0.0015,
                          'bub300 Stainless Steel (for debugging)':0.025}
        pipe_type = st.selectbox(label='Pipe Inner Material',options=list(pipe_materials.keys()))
        pipe_roughness = pipe_materials[pipe_type]/1000 #mm to meters
        
        pipe_diameter = st.number_input(label='Pipe Diameter (in)', value=3., format='%1.2f')
        pipe_diameter = convert.in_to_m(pipe_diameter)
        
        segment_length = st.number_input(label='Segment Length (ft)', value=100)
        segment_length = convert.ft_to_m(segment_length)
        
        number_of_orifices = st.number_input(label='Number of Orifices', value=10, format='%d')
        
        orifice_diameter = st.number_input(label='Orifice Diameter (in)', value=5/8, format='%1.4f')
        orifice_diameter = convert.in_to_m(orifice_diameter)
        
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
    
    
    results = downstream_solve(manifold,
                           air_pressure, 
                           water_pressure, 
                           starting_mdot, 
                           air_temp,
                           verbose=False)
    
    mdot = sum([orifice.mdot for orifice in results.orifices])
    
    st.write('Total Mass flow rate: {:.2f} kg/s'.format(mdot))
    
    atm_pressure = convert.pressure_to_absolute(0)
    rho_air = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
    
    st.write('Total Flow Rate (SI): {:.2f} SCMS @1atm and 68F'.format(air.Q(mdot,rho_air)))
    
    st.write('Total Flow Rate (IMP): {:,.2f} SCFM @1atm and 68F'.format(convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot,rho_air)))))
    
    rho_air = air.rho_air(atm_pressure, T=0) #0C is Standard Temperature in SCFM in BUB300...
    st.write('Total Flow Rate (IMP): {:,.2f} SCFM @1atm and 32F'.format(convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot,rho_air)))))
    
    st.write('Mass Flow Per Orifice (kg/s)')
    st.write([o.mdot for o in results.orifices])
    st.write('Pipe Segment Velocities (m/s)')
    st.write([p.velocity(p.Q) for p in results.pipes])
    st.write('Pipe Segment Pressure Drops (psi)')
    st.write([convert.Pa_to_psi(p.pressure_drop) for p in results.pipes])
    st.write('Total Pipe Pressure Drop (psi)')
    total_p_drop = convert.Pa_to_psi(sum([p.pressure_drop for p in results.pipes]))
    st.write(total_p_drop)
    st.write('Ratio of Total Pressure Drop to Input Pressure')
    st.write(total_p_drop/convert.Pa_to_psi(convert.pressure_to_gauge(air_pressure)))
    st.write('Coefficient of Uniformity')
    st.write(results.orifices[-1].mdot/results.orifices[0].mdot)
    st.write('Pipe Segment Mach Numbers')
    st.write([p.mach_number for p in results.pipes])
    st.write('Pipe Segment Reynolds Number')
    st.write([p.Re for p in results.pipes])
    
    return manifold