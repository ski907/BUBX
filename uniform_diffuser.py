import streamlit as st

from physics import air
from conversions import convert

from geometry.segments import generate_equal_spaced_uniform
from geometry.segments import add_supply_line
from geometry.segments import set_orifice_elevations_on_profile

from solver import downstream_solve
from solver import specify_airflow_get_pressure
from geometry.parse_results import parse_results

def uniform_diffuser():
    

    
    col1, col2, col3, col4 = st.beta_columns(4)
    
    col2.write('Diffuser Geometry')
    
    o_spacing_method = col2.radio('Specify Orifice Placement Method', ('Number','Spacing'))
    
    with col2.form(key='system_geometry_form'):
    
        st.write('Diffuser Geometry Input:')
        #pipe_roughness = st.number_input(label='Pipe Roughness Height (m)', value=0.000025, format='%1.10f')
        
        pipe_materials = {'stainless steel (0.015mm) ':0.015,
                          'galvinized steel (0.15mm)':0.15,
                          'PVC: (0.0015mm)':0.0015,
                          'bub300 Stainless Steel: (0.025mm)':0.025,
                          'Wire Reinforced Rubber Hose: (0.3mm)':0.3}
        pipe_type = st.selectbox(label='Pipe Inner Material (Absolute Roughness)',options=list(pipe_materials.keys()))
        pipe_roughness = pipe_materials[pipe_type]/1000 #mm to meters
        
        pipe_diameter = st.number_input(label='Pipe Diameter (in)', value=3., format='%1.2f')
        pipe_diameter = convert.in_to_m(pipe_diameter)
        
        segment_length = st.number_input(label='Segment Length (ft)', value=100)
        segment_length = convert.ft_to_m(segment_length)
        

        
        if o_spacing_method == 'Number':
            number_of_orifices = st.number_input(label='Number of Orifices', value=10, format='%d')
            spacing = convert.m_to_ft(segment_length/number_of_orifices)
            st.write('Orifice Spacing = {:.2f} ft'.format(spacing))
        
        if o_spacing_method == 'Spacing':
            spacing = st.number_input(label='Orifice Spacing (ft)', value=10., format='%1.2f')
            number_of_orifices = int(segment_length/convert.ft_to_m(spacing))+1
            st.write('Number of Orifices = {}'.format(number_of_orifices))
        
        orifice_diameter = st.number_input(label='Orifice Diameter (in)', value=5/8, format='%1.4f')
        orifice_diameter = convert.in_to_m(orifice_diameter)
        
        sg_submit_button = st.form_submit_button(label='Update/Run')
   
    col1.write('Supply Line Geometry')     
    with col1.form(key='supply_geometry_form'):
        st.write('Supply Line Geometry Input:')
        supply_pipe_type = st.selectbox(label='Supply Pipe Inner Material',options=list(pipe_materials.keys()))
        supply_pipe_roughness = pipe_materials[supply_pipe_type]/1000 #mm to meters
    
        supply_pipe_diameter = st.number_input(label='Supply Pipe Diameter (in)', value=3., format='%1.2f')
        supply_pipe_diameter = convert.in_to_m(supply_pipe_diameter)
    
        supply_pipe_length = st.number_input(label='Supply Pipe Length (ft)', value=100)
        supply_pipe_length = convert.ft_to_m(supply_pipe_length)
        
        supply_submit_button = st.form_submit_button(label='Update/Run')
    
    with col1.form(key='diffuser_elevation_form'):
        st.write("Diffuser Elevation Profile")
        profile = st.text_input("Elevation Profile [offset1,elev1],[offset2,elev2],...", f'[0,0],[{convert.m_to_ft(segment_length)},0]')
        profile = list(eval(profile))

        datum = st.number_input('Datum (ft)',0)

        elevation_submit_button = st.form_submit_button(label='Update/Run')
        
        for point in profile:
            point[0]=convert.ft_to_m(point[0])
            point[1]=convert.ft_to_m(point[1])
            
        datum = convert.ft_to_m(datum)
        
        
    system_geom = generate_equal_spaced_uniform(segment_length, 
                                  pipe_roughness,
                                  pipe_diameter,
                                  number_of_orifices,
                                  orifice_diameter
                                  )
        
    system_geom = set_orifice_elevations_on_profile(system_geom, profile, datum)
        
    system_geom =  add_supply_line(system_geom,
                                    supply_pipe_length,
                                    supply_pipe_roughness,
                                    supply_pipe_diameter)
  
        

    
    col3.write('Boundary Conditions')
    bc_method = col3.radio('Specify Boundary Condition Method', ('Pressure Specified','Airflow Specified'))
    
    starting_mdot = 0.0001
    
    if bc_method == 'Pressure Specified':
    
        air_pressure, water_pressure, air_temp = boundary_conditions_specify_pressure(col3)
        
    

    
    if bc_method == 'Airflow Specified':
        airflow_SCMM, water_pressure, air_temp = boundary_conditions_specify_airflow(col3)
        
        air_pressure = specify_airflow_get_pressure(system_geom, airflow_SCMM, water_pressure, starting_mdot, air_temp)
        
        col3.write(convert.Pa_to_psi(convert.pressure_to_gauge(air_pressure)))


    solved_geom = downstream_solve(system_geom,
                   air_pressure, 
                   water_pressure, 
                   starting_mdot, 
                   air_temp,
                   verbose=False)
 

        
    
    col4.write('Results:')
    mdot = parse_results(solved_geom).total_mdot()
    
    atm_pressure = convert.pressure_to_absolute(0)
    rho_air = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
    
    col4.write('Total Flow Rate (SI): {:.2f} SCMM @1atm and 20C'.format(convert.CMS_to_CMM(air.Q(mdot,rho_air))))
    
    col4.write('Total Flow Rate (IMP): {:,.2f} SCFM @1atm and 68F'.format(convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot,rho_air)))))
    
    col4.write('Air Pressure: {:,.2f} psi'.format(convert.Pa_to_psi(convert.pressure_to_gauge(air_pressure))))

    col4.write('Airflow Per Unit Length: {:,.2f} SCMM/m'.format(parse_results(solved_geom).airflow_per_unit_length()))
    
    col4.write('Estimated Horizontal Surface Velocity (Haehnel 2016): {:,.2f} cm/s'.format(parse_results(solved_geom).horizontal_surface_vel_haehnel2016()))
    col4.write('Estimated Horizontal Surface Velocity (Haehnel 2016): {:,.2f} ft/s'.format(convert.m_to_ft(parse_results(solved_geom).horizontal_surface_vel_haehnel2016())))
    
    coeff_uniformity = parse_results(solved_geom).coefficient_of_uniformity()
    col4.write('Coefficient of Uniformity: {:,.2f}'.format(coeff_uniformity))
    
    if coeff_uniformity < 0.9:
        col4.write('Warning: Coefficient of Uniformity should be greater than 0.9. Imbalance in air screen perfomance likely')
    
    orifice_diffuser_ratio = parse_results(solved_geom).orifice_to_diffuser_area_ratio()
    col4.write('Ratio of Total Orifice Area to Diffuser Area: {:,.2f}'.format(orifice_diffuser_ratio))
    
    if orifice_diffuser_ratio > 0.25:
        col4.write('Warning: Ratio Execeeds Recommended Value of 0.25. Increase Diffuser Size, or Decrease Number/Size of Orifices')
    
    show_detailed_output = col4.checkbox(label='Show Detailed Output?')
    
    if show_detailed_output:
        
        col4.write('Total Mass flow rate: {:.2f} kg/s'.format(mdot))
        
        col4.write('Total Pipe Pressure Drop (psi)')
        total_p_drop = convert.Pa_to_psi(parse_results(solved_geom).pipe_total_pressure_drop())
        col4.write(total_p_drop)
        
        col4.write('Ratio of Total Pressure Drop to Input Pressure')
        col4.write(parse_results(solved_geom).ratio_pressure_drop_to_input())
        
        rho_air = air.rho_air(atm_pressure, T=0) #0C is Standard Temperature in SCFM in BUB300...
        col4.write('Total Flow Rate (IMP): {:,.2f} SCFM @1atm and 32F (temp used by BUB300)'.format(convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot,rho_air)))))
        
        col4.write('Mass Flow Per Orifice (kg/s)')
        col4.write(parse_results(solved_geom).orifice_mdots())
        
        col4.write('Pipe Segment Velocities (m/s)')
        col4.write(parse_results(solved_geom).pipe_velocities())
        
        col4.write('Pipe Segment Pressure Drops (psi)')
        col4.write(list(map(convert.Pa_to_psi,parse_results(solved_geom).pipe_pressure_drops())))
        
        col4.write('Pipe Segment Mach Numbers')
        col4.write(parse_results(solved_geom).pipe_mach_numbers())
        
        col4.write('Pipe Segment Reynolds Number')
        col4.write(parse_results(solved_geom).pipe_reynolds_numbers())
    
    return system_geom

def boundary_conditions_specify_pressure(col):
    with col.form(key='Boundary Conditions'):
        air_pressure = st.number_input(label='Air Pressure (psi)', value=90., format='%1.2f')
        air_pressure = convert.psi_to_Pa(air_pressure)
        air_pressure = convert.pressure_to_absolute(air_pressure)
        
        water_depth = st.number_input(label='Water Depth (ft)', value=30.0, format='%1.1f')
        water_depth = convert.ft_to_m(water_depth)
        water_pressure = convert.pressure_to_absolute(convert.H_m_to_Pa(water_depth))
        
        air_temp = st.number_input(label='Air Temp (C)', value=0, format='%d')
        
        #starting_mdot = st.number_input(label='Starting Airflow (kg/s)', value=0.01)
        
        bc_submit_button = st.form_submit_button(label='Update/Run')
    
    # st.form_submit_button returns True upon form submit
    if bc_submit_button:
        col.write(f'Boundary Conditions Updated')
        
    return air_pressure, water_pressure, air_temp

def boundary_conditions_specify_airflow(col):
    with col.form(key='Boundary Conditions'):
        airflow = st.number_input(label='Air Flow Rate (SCFM)', value=1600, format='%d')
        airflow_SCMM = convert.CMS_to_CMM(convert.CFS_to_CMS(convert.CFM_to_CFS(airflow)))
        
        water_depth = st.number_input(label='Water Depth (ft)', value=30.0, format='%1.1f')
        water_depth = convert.ft_to_m(water_depth)
        water_pressure = convert.pressure_to_absolute(convert.H_m_to_Pa(water_depth))
        
        air_temp = st.number_input(label='Air Temp (C)', value=0, format='%d')
        
        #starting_mdot = st.number_input(label='Starting Airflow (kg/s)', value=0.01)
        
        bc_submit_button = st.form_submit_button(label='Update/Run')
    
    # st.form_submit_button returns True upon form submit
    if bc_submit_button:
        col.write(f'Boundary Conditions Updated')
        
    return airflow_SCMM, water_pressure, air_temp