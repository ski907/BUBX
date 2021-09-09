from geometry import Pipe
from geometry import Orifice

class Segment:
    def __init__(self, features):
        self.features = features
        
        self.orifices = [feature for feature in features if isinstance(feature, Orifice)]
        self.pipes = [feature for feature in features if isinstance(feature, Pipe)]
        