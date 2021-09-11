import streamlit as st

import numpy as np

#from scipy.optimize import minimize


from conversions import convert


from uniform_diffuser import *
import plots

st.set_page_config(layout="wide")
st.title('BUBX version 0.0')
st.write("An application for solving air demands for submerged manifolds")

system_geometry = uniform_diffuser()

col4, col5 = st.beta_columns((1,2))

with col4.form(key='Flow Range Plot Parameters'):
    
    min_depth = st.number_input('Minimum Water Depth (ft)',value = 30.)
    max_depth = st.number_input('Maximum Water Depth (ft)',value = 50.)
    
    
    
    min_ap = st.number_input('Minimum Air Pressure (psi)',value = 25)
    max_ap = st.number_input('Maximum Air Pressure (psi)',value = 100)
    air_press_inc = st.number_input('Air Pressure Increments (psi)',value = 10)
    
    units = st.radio('Units',['SI','English'])
    
    air_pressures = range(int(min_ap),int(max_ap+1),int(air_press_inc))
    
    #convert everything to SI
    min_depth_m = convert.ft_to_m(min_depth) 
    max_depth_m = convert.ft_to_m(max_depth)      
    depth_lims_m = (min_depth_m,max_depth_m)
    air_pressures_Pa = [convert.pressure_to_absolute(convert.psi_to_Pa(pressure)) for pressure in air_pressures]
    
    fr_submit_button = st.form_submit_button(label='Update')
    
    # st.form_submit_button returns True upon form submit
    if fr_submit_button:
        
        col5.write(plots.plot_flow_vs_depth(system_geometry,
                       depth_lims_m,
                       air_pressures_Pa,
                       air_temp=0,
                       units=units))


