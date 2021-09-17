import sys
import math
import numpy as np
#from scipy.optimize import minimize

from physics import air
from conversions import convert

from geometry import Pipe
from geometry import Orifice

from geometry.segments import generate_equal_spaced_uniform
from geometry.segments import generate_point
from geometry.segments import add_supply_line
from geometry.segments import generate_repeating_segments
from geometry.parse_results import parse_results


#def mdot_calc(starting_mdot, air_pressure, water_pressure, system_geometry, air_temp, mode='solve'):    
#    '''Function to be optimized with scipy solvers 
#    given mdot input, returns residual mdot at downstream end
#    The system is solved when the residual mdot returned in close to zero'''
#    
#    mdot_calc = abs(starting_mdot)
#        
#    air_pressure_calc = air_pressure
#
#    
#    for feature in system_geometry.features:
#        
#        if isinstance(feature, Orifice):
#            
#            feature.upstream_pressure = air_pressure_calc
#            feature.downstream_pressure = water_pressure
#                        
#            mdot_calc -= feature.mdot
#            
#        
#        if isinstance(feature, Pipe):
#            
#            feature.rho_air = air.rho_air(air_pressure_calc, air_temp)
#            feature.Q = air.Q(mdot_calc, feature.rho_air)
#                        
#            air_pressure_calc -= feature.pressure_drop
#            
#            if air_pressure_calc < 0:
#                mdot_calc = 1000
#    
#    if mode == 'solve':        
#        return abs(mdot_calc)
#    elif mode == 'return system geometry':
#        return system_geometry
#    
#    
#def scipy_minimize_solve(segment, air_pressure, water_pressure, starting_mdot, air_temp):
#    '''This approach is not robust yet- in test cases it converges on 
#    incorrect values i.e. close to zero, or very large values'''
#    
#    
#    scipy_results = minimize(mdot_calc,
#                             starting_mdot,
#                             args=(air_pressure, water_pressure, segment, air_temp),
#                             method='Nelder-Mead', 
#                             tol=1e-4)
#    
#    #print(scipy_results)
#    print('scipy-minimize Converged on {} kg/s after {} iterations'.format(abs(scipy_results.x[0]), scipy_results.nit))
#    
#    segment = mdot_calc(abs(scipy_results.x[0]), 
#                         air_pressure, 
#                         water_pressure,
#                         segment,
#                         air_temp, 
#                         mode='return system geometry')
#    
#    return segment
    

def are_within_tolerance(a, b, tolerance=1e-10):
    return math.fabs(a-b) < tolerance

def downstream_solve(system_geometry, air_pressure, water_pressure, starting_mdot, air_temp, verbose=False):
    #perhaps update this to bisection solver. It breaks for air pressures near the water pressure
    if water_pressure > air_pressure:
        sys.exit('Water Pressure Exceeds Internal Air Pressure')
    
    mdot_calc = starting_mdot
    iterations = 0 #initialize
    converged = False
    run_broken = True
    max_iterations = 1000
    mdot_increase = 10 #percent
    mdot_decrease = 10 #percent
    
    
    while not converged or run_broken:
        
        mdot_calc = starting_mdot
        
        air_pressure_calc = air_pressure
        iterations += 1
    
        for feature in system_geometry.features:
            
            if isinstance(feature, Orifice):
                
                feature.temp = air_temp
                feature.upstream_pressure = air_pressure_calc
                feature.downstream_pressure = water_pressure
                            
                mdot_calc -= feature.mdot
                
                if mdot_calc < 0:
                    if verbose: print('Orifice mass rate out exceeds available air mass rate, increasing air flow {}%'.format(mdot_increase))
                    starting_mdot *= (1 + mdot_increase/100)
                    run_broken = True
                    break
            
            if isinstance(feature, Pipe):
                feature.air_temp = air_temp
                method = "average pressure"
                #method = "upstream pressure"
                
                
                if method == "average pressure":
                    feature.p1 = air_pressure_calc
                    feature.mdot = mdot_calc
                    air_pressure_calc -= feature.pressure_drop_average

                if method == "upstream pressure":               
                    feature.rho_air = air.rho_air(air_pressure_calc, air_temp)
                    feature.Q = air.Q(mdot_calc, feature.rho_air)
                    air_pressure_calc -= feature.pressure_drop
                
                if air_pressure_calc < water_pressure:
                    if verbose: print('Air pressure drop along pipe results in pressure below external hydrostatic pressure, decreasing air flow by {}% of residual mass flow'.format(mdot_decrease))
                    starting_mdot -= mdot_calc * (mdot_decrease/100)
                    run_broken = True
                    break
            
            run_broken = False
        
        if are_within_tolerance(mdot_calc, 0, tolerance = 0.000001):
            converged = True
        
        else:
            if verbose: print('Airflow mass rate in exceeds mass rate out, decreasing airflow by {}% of residual mass flow'.format(mdot_decrease))
            starting_mdot -= mdot_calc * (mdot_decrease/100)
            if verbose: print('Mass rate = {} kg/s'.format(starting_mdot))
            
        if iterations > max_iterations:
            sys.exit('Max iterations of {} exceeded'.format(max_iterations))
            
        
    if verbose: print('Converged on {} kg/s after {} iterations'.format(starting_mdot, iterations))
    
    atm_pressure = convert.pressure_to_absolute(0)
    rho_air = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
    
    if verbose: print('This is equal to {} CMS'.format(air.Q(starting_mdot,rho_air)))
    if verbose: print('This is equal to {} CFM'.format(convert.CFS_to_CFM(convert.CMS_to_CFS(air.Q(starting_mdot,rho_air)))))
    
    if verbose: print([o.mdot for o in system_geometry.orifices])

    return system_geometry

def specify_airflow_get_pressure_OLD(system_geometry, airflow_SCMM, water_pressure, starting_mdot, air_temp):
    #this works OK, but gets flaky when run at low airflows see new method which uses method of bisection
    atm_pressure = convert.pressure_to_absolute(0)
    rho_air_standard = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
    mdot_target = air.mdot(convert.CMM_to_CMS(airflow_SCMM), rho_air_standard)
    
    air_pressure = water_pressure * 1.01
    
    iterations = 0
    converged = False
    while not converged:
        iterations += 1
        solved_geo  = downstream_solve(system_geometry, air_pressure, water_pressure, starting_mdot, air_temp)
        mdot_calc = parse_results(solved_geo).total_mdot()
        
        if are_within_tolerance(mdot_target, mdot_calc, tolerance = 0.000001):
            converged = True
        else:
            if mdot_calc < mdot_target:
                air_pressure += air_pressure * (mdot_target-mdot_calc)/mdot_target * 1
                
            elif mdot_calc > mdot_target:
                air_pressure -= air_pressure * (mdot_calc-mdot_target)/mdot_target * 0.3
        print(air_pressure)
    print(iterations)            
    return air_pressure

def specify_airflow_get_pressure(system_geometry, airflow_SCMM, water_pressure, starting_mdot, air_temp, tol = 0.0001):
    #updated solver which uses method of bisection
    max_iterations = 1000
    iterations = 0
    atm_pressure = convert.pressure_to_absolute(0)
    rho_air_standard = air.rho_air(atm_pressure, T=convert.F_to_C(68)) #68F is Standard Temperature in SCFM in north america
    mdot_target = air.mdot(convert.CMM_to_CMS(airflow_SCMM), rho_air_standard)
    
    max_pressure = 1000
    
    a = water_pressure * 1.01
    b = convert.psi_to_Pa(max_pressure)
    
    def f(air_pressure):
        solved_geo  = downstream_solve(system_geometry, air_pressure, water_pressure, starting_mdot, air_temp)
        mdot_calc = parse_results(solved_geo).total_mdot()     
        return mdot_calc - mdot_target
    
    f_a = f(a)
    f_b = f(b)
    
    
    if np.sign(f_a) == np.sign(f_b):
        raise Exception(
         'The scalars a and b do not bound a root. The airflow selected either exceeds that which can be supplied by {} psi, or is so low that the associated pressure is close to that of the water pressure'.format(max_pressure))
    

    #basic bisection method solver
    converged = False
    while not converged:    
        iterations += 1
        
        m = (a + b)/2
        f_m = f(m)
        
        print(m)
        if np.abs(f_m) < tol:
            return m
        elif np.sign(f_a) == np.sign(f_m):
            #m is an improvement on a
            a = m
            f_a = f(a)
        elif np.sign(f_b) == np.sign(f_m):
            #m is an improvement on b
            b = m
            f_b = f(b)
        
        if iterations > max_iterations:
            sys.exit('Max iterations of {} exceeded (In specify airflow get pressure routine)'.format(max_iterations))


def main():
    
    manifold = generate_equal_spaced_uniform(segment_length=100, 
                                  pipe_roughness=0.000025,
                                  pipe_diameter= 0.0762009266,
                                  number_of_orifices=100,
                                  orifice_diameter=0.0015875193
                                  )
    
#    manifold = generate_point(segment_length = 9.144,
#                              pipe_roughness = 0.000025,
#                              pipe_diameter = 0.0762009266,
#                              orifice_diameter = 0.0127
#                              )
#   
#    manifold = generate_repeating_segments(segment_length=20, 
#                                  pipe_roughness=0.000025,
#                                  pipe_diameter= 0.0762009266,
#                                  number_of_orifices=20,
#                                  orifice_diameter=0.0015875193,
#                                  number_of_repeats=5,
#                                  connector_length=5.,
#                                  connector_roughness =0.000025,
#                                  connector_diameter = 0.0762009266
#                                  )
    
    manifold = add_supply_line(manifold, 
                               supply_length = 100, 
                               pipe_roughness = 0.000025, 
                               pipe_diameter = 0.0762009266
                               )
    
    
    
    air_pressure = 100 #psi
    air_pressure = convert.psi_to_Pa(air_pressure)
    air_pressure = convert.pressure_to_absolute(air_pressure)
    
    water_depth = 30 #ft
    water_depth = convert.ft_to_m(water_depth)
    water_pressure = convert.pressure_to_absolute(convert.H_m_to_Pa(water_depth))

    air_temp = 0
        
    starting_mdot = 1
    
    air_flow_SCMM =30
    
    ap = specify_airflow_get_pressure(manifold,air_flow_SCMM, water_pressure, starting_mdot, air_temp)
    
    geom = downstream_solve(manifold,
                               air_pressure, 
                               water_pressure, 
                               starting_mdot, 
                               air_temp,
                               verbose=True)
    
    
    
    print(sum([orifice.mdot for orifice in geom.orifices]))
    
#    scipy_results = scipy_minimize_solve(manifold, 
#                                         air_pressure, 
#                                         water_pressure, 
#                                         starting_mdot, 
#                                         air_temp)
#    
#    print(sum([orifice.mdot for orifice in scipy_results.orifices]))
    
    
    
if __name__ == '__main__':
    main()