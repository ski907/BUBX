from segments import make_segments
from segments import generate_equal_spaced_uniform

pipes = [[10, .2, .01],
         [20, .2, .02],
         [40, .4, .01]]

orifices = [[0.025], [0.05], [0.03]]

seg1 = make_segments(pipes, orifices)

seg1.orifices

seg2 = generate_equal_spaced_uniform(segment_length=100, 
                                  pipe_roughness=.02,
                                  pipe_diameter=2,
                                  number_of_orifices=5,
                                  orifice_diameter=0.0625
                                  )


seg2.orifices

