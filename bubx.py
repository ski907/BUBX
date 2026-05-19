import streamlit as st
import numpy as np
from conversions import convert
from uniform_diffuser import *
import plots
from pdf_report import generate_pdf_report
from datetime import datetime

# Page configuration
st.set_page_config(
    page_title="BUBX - Air Demand Calculator",
    page_icon="🫧",
    layout="wide"
)

# Header section
st.title('🫧 BUBX')
st.markdown("**Version 0.2** - An application for solving air demands for submerged manifolds")
st.markdown("---")

# Initialize system geometry
#system_geometry = uniform_diffuser()

# Main content area with tabs
tab_input, tab_results, tab_flow_plots, tab_flow_range_plots, tab_about = st.tabs(["⚙️ Inputs", "🔢 Results", "📊 Orifice Flow Plots", "📈 Flow Range Analysis", "ℹ️ About"])

with tab_input:
    col1, col2, col3 = st.columns(3)
    
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
        
        pipe_diameter_in = st.number_input(label='Pipe Diameter (in)', value=3., format='%1.2f')
        pipe_diameter = convert.in_to_m(pipe_diameter_in)

        segment_length_ft = st.number_input(label='Segment Length (ft)', value=100)
        segment_length = convert.ft_to_m(segment_length_ft)
        

        
        if o_spacing_method == 'Number':
            number_of_orifices = st.number_input(label='Number of Orifices', value=10, format='%d')
            spacing = convert.m_to_ft(segment_length/number_of_orifices)
            st.write('Orifice Spacing = {:.2f} ft'.format(spacing))
        
        if o_spacing_method == 'Spacing':
            spacing = st.number_input(label='Orifice Spacing (ft)', value=10., format='%1.2f')
            number_of_orifices = int(segment_length/convert.ft_to_m(spacing))+1
            st.write('Number of Orifices = {}'.format(number_of_orifices))
        
        orifice_diameter_in = st.number_input(label='Orifice Diameter (in)', value=5/8, format='%1.4f')
        orifice_diameter = convert.in_to_m(orifice_diameter_in)
        
        sg_submit_button = st.form_submit_button(label='Update/Run')
   
    col1.write('Supply Line Geometry')     
    with col1.form(key='supply_geometry_form'):
        st.write('Supply Line Geometry Input:')
        supply_pipe_type = st.selectbox(label='Supply Pipe Inner Material',options=list(pipe_materials.keys()))
        supply_pipe_roughness = pipe_materials[supply_pipe_type]/1000 #mm to meters
    
        supply_pipe_diameter_in = st.number_input(label='Supply Pipe Diameter (in)', value=3., format='%1.2f')
        supply_pipe_diameter = convert.in_to_m(supply_pipe_diameter_in)

        supply_pipe_length_ft = st.number_input(label='Supply Pipe Length (ft)', value=100)
        supply_pipe_length = convert.ft_to_m(supply_pipe_length_ft)
        
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
        
        
    system_geometry = generate_equal_spaced_uniform(segment_length, 
                                  pipe_roughness,
                                  pipe_diameter,
                                  number_of_orifices,
                                  orifice_diameter
                                  )
        
    system_geometry = set_orifice_elevations_on_profile(system_geometry, profile, datum)
        
    system_geometry =  add_supply_line(system_geometry,
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
        
        air_pressure = specify_airflow_get_pressure(system_geometry, airflow_SCMM, water_pressure, starting_mdot, air_temp)


    solved_geom = downstream_solve(system_geometry,
                   air_pressure,
                   water_pressure,
                   starting_mdot,
                   air_temp,
                   verbose=False)

    report_inputs = {
        'pipe_type': pipe_type,
        'pipe_diameter_in': pipe_diameter_in,
        'segment_length_ft': segment_length_ft,
        'number_of_orifices': number_of_orifices,
        'orifice_diameter_in': orifice_diameter_in,
        'supply_pipe_type': supply_pipe_type,
        'supply_pipe_diameter_in': supply_pipe_diameter_in,
        'supply_pipe_length_ft': supply_pipe_length_ft,
        'bc_method': bc_method,
        'water_depth_ft': convert.m_to_ft(convert.Pa_to_H_m(convert.pressure_to_gauge(water_pressure))),
    }

with tab_results:
    st.header('Results')
    
    # Calculate base values once
    results = parse_results(solved_geom)
    mdot = results.total_mdot()
    atm_pressure = convert.pressure_to_absolute(0)
    
    # Standard conditions for different regions
    rho_air_68f = air.rho_air(atm_pressure, T=convert.F_to_C(68))  # 68F - North American standard
    rho_air_0c = air.rho_air(atm_pressure, T=0)  # 0C - BUB300 standard
    
    # Main Results Table
    st.subheader('Flow Rates & System Parameters')
    
    # Calculate all values
    flow_rate_si = convert.CMS_to_CMM(air.Q(mdot, rho_air_68f))
    flow_rate_imp = convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot, rho_air_68f)))
    air_pressure_psi = convert.Pa_to_psi(convert.pressure_to_gauge(air_pressure))
    airflow_per_length = results.airflow_per_unit_length()
    airflow_per_length_imp = convert.CMM_to_CFM(airflow_per_length)
    surface_vel_cms = results.horizontal_surface_vel_haehnel2016()
    surface_vel_fts = convert.m_to_ft(surface_vel_cms)
    
    # Create results table
    import pandas as pd
    results_data = {
        'Parameter': [
            'Total Flow Rate',
            'Air Pressure', 
            'Airflow Per Unit Length',
            'Surface Velocity (Haehnel 2016)'
        ],
        'SI Units': [
            f'{flow_rate_si:.2f} SCMM @ 1atm, 20°C',
            f'{convert.pressure_to_gauge(air_pressure)/1000000:,.2f} MPa',
            f'{airflow_per_length:,.2f} SCMM/m',
            f'{surface_vel_cms:,.2f} m/s'
        ],
        'Imperial Units': [
            f'{flow_rate_imp:,.2f} SCFM @ 1atm, 68°F',
            f'{air_pressure_psi:,.2f} psi',
            f'{airflow_per_length_imp:,.2f} SCFM/m', 
            f'{surface_vel_fts:,.2f} ft/s'
        ]
    }
    
    df_results = pd.DataFrame(results_data)
    st.dataframe(df_results, use_container_width=True, hide_index=True)
    
    # Performance Metrics Table
    st.subheader('Performance Metrics')
    
    coeff_uniformity = results.coefficient_of_uniformity()
    orifice_diffuser_ratio = results.orifice_to_diffuser_area_ratio()
    
    # Create performance table
    performance_data = {
        'Metric': [
            'Coefficient of Uniformity',
            'Orifice/Diffuser Area Ratio'
        ],
        'Value': [
            f'{coeff_uniformity:.3f}',
            f'{orifice_diffuser_ratio:.3f}'
        ],
        'Target': [
            '> 0.9',
            '≤ 0.25'
        ],
        'Status': [
            '✅ OK' if coeff_uniformity >= 0.9 else '⚠️ Warning',
            '✅ OK' if orifice_diffuser_ratio <= 0.25 else '⚠️ Warning'
        ]
    }
    
    df_performance = pd.DataFrame(performance_data)
    st.dataframe(df_performance, use_container_width=True, hide_index=True)
    
    # Show warnings below table
    if coeff_uniformity < 0.9:
        st.warning("Coefficient of Uniformity should be > 0.9. Air screen imbalance likely.")
    
    if orifice_diffuser_ratio > 0.25:
        st.warning("Ratio exceeds recommended 0.25. Consider increasing diffuser size or reducing orifice count/size.")

    # PDF Export
    st.divider()
    if st.button("Generate PDF Report", type="primary"):
        with st.spinner("Building PDF..."):
            pdf_bytes = generate_pdf_report(
                solved_geom, report_inputs, air_pressure, water_pressure, air_temp
            )
        st.download_button(
            label="Download PDF",
            data=pdf_bytes,
            file_name=f"BUBX_Report_{datetime.now().strftime('%Y%m%d_%H%M')}.pdf",
            mime="application/pdf",
        )

    # Detailed Output Toggle
    st.divider()
    show_detailed_output = st.checkbox("Show Detailed Output", help="Display additional technical parameters")
    
    if show_detailed_output:
        st.subheader('Detailed Technical Output')
        
        # Flow and Pressure Details Table
        with st.expander("Flow and Pressure Analysis", expanded=True):
            
            total_p_drop_psi = convert.Pa_to_psi(results.pipe_total_pressure_drop())
            pressure_ratio = results.ratio_pressure_drop_to_input()
            flow_rate_bub300 = convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(mdot, rho_air_0c)))
            
            detailed_data = {
                'Parameter': [
                    'Total Mass Flow Rate',
                    'Total Pipe Pressure Drop', 
                    'Pressure Drop Ratio',
                    'Flow Rate (BUB300 Standard)'
                ],
                'Value': [
                    f'{mdot:.2f} kg/s',
                    f'{total_p_drop_psi:.2f} psi',
                    f'{pressure_ratio:.3f}',
                    f'{flow_rate_bub300:,.2f} SCFM'
                ],
                'Notes': [
                    '',
                    'Total pressure loss through pipe system',
                    'Total drop / Input pressure',
                    '@ 1atm and 32°F (BUB300 temp standard)'
                ]
            }
            
            df_detailed = pd.DataFrame(detailed_data)
            st.dataframe(df_detailed, use_container_width=True, hide_index=True)
        
        # Per-Orifice Analysis
        with st.expander("Per-Orifice Analysis"):
            st.write("**Mass Flow Per Orifice (kg/s)**")
            orifice_flows = results.orifice_mdots()
            st.dataframe(orifice_flows, use_container_width=True)
        
        # Pipe Segment Analysis
        with st.expander("Pipe Segment Analysis"):
            segment_data = {
                'Velocities (m/s)': results.pipe_velocities(),
                'Pressure Drops (psi)': [convert.Pa_to_psi(p) for p in results.pipe_pressure_drops()],
                'Mach Numbers': results.pipe_mach_numbers(),
                'Reynolds Numbers': results.pipe_reynolds_numbers()
            }
            
            # Create a more readable table format
            import pandas as pd
            df = pd.DataFrame(segment_data)
            df.index.name = 'Segment'
            st.dataframe(df, use_container_width=True)


with tab_flow_plots:
    st.header("Orifice Flow Analysis")
    
    # Unit selection in a more compact way
    col_unit1, col_unit2, col_unit3 = st.columns([1, 1, 3])
    with col_unit1:
        orifice_plot_units = st.radio(
            'Display Units:', 
            ('SI', 'English'),
            horizontal=True,
            key='orifice_units'
        )
    
    # Display plots with better spacing
    st.subheader("Orifice Flow Distribution")
    st.write(plots.plot_orifice_flows(system_geometry, units=orifice_plot_units))
    
    st.subheader("Horizontal Velocity Profile")
    st.write(plots.plot_horizontal_velocities(system_geometry, units=orifice_plot_units))

with tab_flow_range_plots:
    st.header("Flow vs. Depth Analysis")
    
    # Create two columns for better layout
    col_params, col_plot = st.columns([1, 2])
    
    with col_params:
        st.subheader("Analysis Parameters")
        
        with st.form(key='flow_range_params'):
            # Water Depth Section
            st.markdown("#### Water Depth Range")
            col_depth1, col_depth2 = st.columns(2)
            with col_depth1:
                min_depth = st.number_input(
                    'Minimum (ft)',
                    value=30.0,
                    min_value=0.0,
                    step=5.0,
                    help="Minimum water depth for analysis"
                )
            with col_depth2:
                max_depth = st.number_input(
                    'Maximum (ft)',
                    value=50.0,
                    min_value=0.0,
                    step=5.0,
                    help="Maximum water depth for analysis"
                )
            
            # Air Pressure Section
            st.markdown("#### Air Pressure Range")
            min_ap = st.number_input(
                'Minimum Pressure (psi)',
                value=25,
                min_value=0,
                step=5,
                help="Starting air pressure"
            )
            max_ap = st.number_input(
                'Maximum Pressure (psi)',
                value=100,
                min_value=0,
                step=5,
                help="Ending air pressure"
            )
            air_press_inc = st.number_input(
                'Pressure Increment (psi)',
                value=10,
                min_value=1,
                step=5,
                help="Step size for pressure range"
            )
            
            # Display Units
            st.markdown("#### Display Options")
            units = st.radio(
                'Output Units',
                ['SI', 'English'],
                horizontal=True,
                help="Units for the generated plot"
            )
            
            # Submit button with styling
            submit_button = st.form_submit_button(
                label='🔄 Update Analysis',
                use_container_width=True,
                type='primary'
            )
    
    with col_plot:
        if submit_button or 'flow_plot_generated' not in st.session_state:
            # Calculate air pressure range
            air_pressures = range(int(min_ap), int(max_ap + 1), int(air_press_inc))
            
            # Convert to SI units for calculation
            min_depth_m = convert.ft_to_m(min_depth)
            max_depth_m = convert.ft_to_m(max_depth)
            depth_lims_m = (min_depth_m, max_depth_m)
            air_pressures_Pa = [
                convert.pressure_to_absolute(convert.psi_to_Pa(pressure)) 
                for pressure in air_pressures
            ]
            
            # Display plot
            st.subheader("Flow vs. Depth Results")
            
            # Show current parameters
            with st.expander("Current Analysis Parameters", expanded=True):
                param_col1, param_col2 = st.columns(2)
                with param_col1:
                    st.metric("Depth Range", f"{min_depth:.1f} - {max_depth:.1f} ft")
                    st.metric("Pressure Range", f"{min_ap} - {max_ap} psi")
                with param_col2:
                    st.metric("Pressure Steps", f"{air_press_inc} psi")
                    st.metric("Number of Curves", len(air_pressures))
            
            # Generate and display plot
            st.write(
                plots.plot_flow_vs_depth(
                    system_geometry,
                    depth_lims_m,
                    air_pressures_Pa,
                    air_temp=0,
                    units=units
                )
            )
            
            st.session_state['flow_plot_generated'] = True

with tab_about:
    st.header("About BUBX")
    
    col_about1, col_about2 = st.columns(2)
    
    with col_about1:
        st.markdown("""
        ### Application Overview
        BUBX is a specialized tool for analyzing air demand in submerged manifold systems.
        
        **Key Features:**
        - Orifice flow distribution analysis
        - Horizontal velocity profiling
        - Flow vs. depth relationship modeling
        - Support for both SI and English units
        
        ### How to Use
        1. **Orifice Analysis Tab**: View the orifice flow distribution and velocity profiles
        2. **Flow Range Analysis Tab**: Configure parameters to analyze flow behavior across different depths and pressures
        3. Select your preferred units (SI or English) for each analysis
        """)
    
    with col_about2:
        st.markdown("""
        ### Technical Information
        
        **System Components:**
        - Uniform diffuser geometry
        - Pressure-depth relationships
        - Flow rate calculations
        
        **Input Parameters:**
        - Water depth range (ft)
        - Air pressure range (psi)
        - Temperature conditions
        
        ### Version History
        - **v0.1** - Initial release with basic functionality
        - **v0.2** - Updated Interface
        """)
        
        # Add system status or additional info if needed
        with st.expander("System Information"):
            st.info("Current system: Uniform Diffuser")
            st.code(f"Geometry Type: {type(system_geometry).__name__}")

# Footer
st.markdown("---")
st.markdown(
    """
    <div style='text-align: center; color: gray; font-size: 12px;'>
    BUBX v0.2 | Air Demand Calculator for Submerged Manifolds
    </div>
    """,
    unsafe_allow_html=True
)