import math
from physics import air

        
class Orifice:
    
    def __init__(self, diameter: float, contraction_coeff = 0.9):
        self.diameter = diameter
        self.contraction_coeff = contraction_coeff
        self.area = math.pi * self.diameter**2 * (1/4)
        
        self.gamma = 1.4
        self.rc = (2/(self.gamma + 1))**(self.gamma/(self.gamma-1))
        self.temp = 0

        self.upstream_pressure = 0
        self.downstream_pressure = 0
    
    @property        
    def mdot(self):
        
        rho = air.rho_air(self.upstream_pressure, self.temp)
        self.r = self.downstream_pressure / self.upstream_pressure
        
        if self.r > self.rc:
            kn = self.kn_subsonic
            self.choked = False
        else:
            kn = self.kn_sonic
            self.choked = True
        
        return self.contraction_coeff * kn * self.area * math.sqrt(self.upstream_pressure * rho) 
        
    @property
    def kn_subsonic(self):
        return math.sqrt(2*self.gamma/(self.gamma-1) * self.r**(2/self.gamma) \
                         * (1 - self.r ** ((self.gamma - 1)/self.gamma)))
    @property
    def kn_sonic(self):
        return math.sqrt(self.gamma * (2/(self.gamma+1))**((self.gamma+1)/(self.gamma-1)))
    

    
    def Q_at_water_pressure(self, downstream_pressure):
        return air.Q(self.mdot, downstream_pressure)
    


    