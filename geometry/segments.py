from geometry import Pipe
from geometry import Orifice
from geometry import Segment

import operator
from functools import reduce
import itertools
    
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



