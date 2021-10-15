from geometry import Pipe
from geometry import Orifice
from geometry import Segment

import operator
from functools import reduce
import itertools

from scipy import interpolate
    
def generate_equal_spaced_uniform(segment_length: float, 
                                  pipe_roughness: float,
                                  pipe_diameter: float,
                                  number_of_orifices: int,
                                  orifice_diameter: float
                                  ):
    
    number_of_pipes = number_of_orifices-1
    pipe_length = segment_length / (number_of_pipes)
    
    pipes = [Pipe(pipe_length, pipe_diameter, pipe_roughness) for pipe_number in range(number_of_pipes)]
    
    orifices = [Orifice(orifice_diameter) for orifice_number in range(number_of_orifices)]
    
    features = reduce(operator.add, itertools.zip_longest(orifices, pipes))
    
    #this gets rid on None objects if orifice and pipe list lengths are unequal
    features = [feature for feature in features if feature] 
    
    return Segment(features)
    
def generate_irregular_segment(pipe_lengths: list,
                               pipe_roughnesses: list,
                               pipe_diameters: list,
                               orifice_diameters: list
                               ):

    pipes = [Pipe(pipe_length, pipe_diameter, pipe_roughness) 
             for pipe_length, pipe_diameter, pipe_roughness 
             in zip(pipe_lengths, pipe_roughnesses,pipe_diameters)]
    
    orifices = [Orifice(orifice_diameter) for orifice_diameter in orifice_diameters]
    
    features = reduce(operator.add, itertools.zip_longest(orifices, pipes))
    
    #this gets rid on None objects if orifice and pipe list lengths are unequal
    features = [feature for feature in features if feature] 
    
    return Segment(features)

def generate_connector_segment(segment_length: float,
                            pipe_roughness: float,
                            pipe_diameter: float
                            ):    
    
    features = [Pipe(segment_length, pipe_diameter, pipe_roughness)]
    
    return Segment(features)
    
    
def generate_point(segment_length: float,
                   pipe_roughness: float,
                   pipe_diameter: float,
                   orifice_diameter: float
                   ):
    
    features = [Pipe(segment_length, pipe_diameter, pipe_roughness), 
                Orifice(orifice_diameter)]

    return Segment(features)

def generate_repeating_segments(segment_length: float, 
                                pipe_roughness: float,
                                pipe_diameter: float,
                                number_of_orifices: int,
                                orifice_diameter: float,
                                number_of_repeats: int,
                                connector_length: float,
                                connector_roughness: float,
                                connector_diameter: float):
    
    features = []
    
    for __ in range(number_of_repeats):
        
        features.extend(generate_equal_spaced_uniform(segment_length,
                                      pipe_roughness,
                                      pipe_diameter,
                                      number_of_orifices,
                                      orifice_diameter
                                      ).features
                                      )
        
        features.extend(generate_connector_segment(connector_length,
                                               connector_roughness,
                                               connector_diameter
                                               ).features
                                               )
        
    return Segment(features)   

def add_supply_line(segment: Segment,
                   supply_length: float,
                   pipe_roughness: float,
                   pipe_diameter: float,
                   ):
    
    features = segment.features
    features.insert(0,Pipe(supply_length, pipe_diameter, pipe_roughness))
    
    return Segment(features)
    

def set_orifice_elevations_on_profile(segment: Segment, profile: list, datum):
    
    elevations = []
    offsets = []
    
    for point in profile:
        offsets.append(point[0])
        elevations.append(point[1]-datum)
    
    offsets.append(offsets[-1] + 0.1) #need to add a little bump to deal with rounding errors
    elevations.append(elevations[-1])
    
    f = interpolate.interp1d(offsets, elevations)
    
    diffuser_offset = 0
    
    for feature in segment.features:
        if isinstance(feature, Orifice):
            feature.elevation_above_datum = float(f(diffuser_offset))
            feature.datum = datum
        
        if isinstance(feature, Pipe):
            diffuser_offset += feature.length
            feature.elevation_above_datum = float(f(diffuser_offset))
            feature.datum = datum
            
            
    features = segment.features        
    return Segment(features)


        
    
        
    
