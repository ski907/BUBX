import math
import fluids
from physics import air

class Pipe:
    def __init__(self, length: float, diameter: float, epsilon: float):
        self.length = length
        self.diameter = diameter
        self.epsilon = epsilon
        self.area = math.pi * self.diameter**2 * (1/4)
        self.Q = 0
        self.rho_air = 0
        self.p1 = 0
        self.p2 = 0
        self.p_average = 0
        self.F_METHOD = 'clamond'
        self.air_temp = 0
        self.mdot = 0
        #self.F_METHOD = 'BUB300 Legacy'
    
  
    def velocity(self, Q: float):
        return Q/self.area
    
    @property
    def Re(self):
        mu = 17.15*10**-6 #N s/m2 dynamic viscosity of air @ 0C
        v = mu / self.rho_air #kinematic viscosity as a function density (function of T and P), kinematic viscosity weak function of lower pressures
        #v = 1.338*10**-5 #m^2/s kinematic viscosity of air @ 0C @ at 1 atm
        Re = self.velocity(self.Q) * self.diameter / v
        return Re
    
    @property
    def f(self):
        '''Clamond method from fluids library is a lot faster than the 
        iterative solution for f found in BUB300. 
        In testing the results aren't that different
        '''
        if self.F_METHOD == 'clamond':
            self._f = fluids.friction.friction_factor(Re=self.Re, \
                        eD=self.epsilon/self.diameter,Method='Clamond')
        
        if self.F_METHOD == 'BUB300 Legacy':
            self._f = self.f_BUB300_legacy
            
        return self._f  

    
    @property
    def pressure_drop(self): 
        '''right now this is just Darcy-Weisbach
        According to Crane 1999 (TP-410 Flow of Fluids Through Valves, Fittings adn Pipe)
        D-W is OK for compressible flow when pressure drops are less than 40% 
        across a pipe segment which is usually the case for air screens.
        Legacy BUB300 used some sort of correction to the D-W pressure drop which
        looks sort of like an assumption of isentropic conditions (no temp change)
        but it is not clear what the source of the correction equation is.
        '''
        Pdel = self.f*self.rho_air*self.length\
        *self.velocity(self.Q)**2/(2*self.diameter)
        
        #bub300 legacy method for correcting pdrop for compressible flow origin unknown
        #p2 = self.p1 - Pdel
        #Pdel = Pdel * 2 * (1 + ((2*self.diameter)/(self.length * self.f)\
        #    *(math.log(self.p1/p2) / (1 + p2/self.p1)) ))
        
        return Pdel
    
    @property
    def mach_number(self):
        average_pressure = (2 * self.p1-self.pressure_drop_average)/2
        average_density = air.rho_air(average_pressure, self.air_temp)
        speed_of_sound = math.sqrt(1.4 * average_pressure / average_density)
        return self.velocity(self.Q) / speed_of_sound
    
    @property
    def pressure_drop_average(self):
        '''This was updated after reviewing Crane 1999 which suggested
        using the average values of pressure and density in cases where 
        the pressure drop could exceed 10% of the upstream pressure
        vs assuming the upstream pressure is appropriate to represent the whole pipe.
        This method assumes that the mass flow rate mdot and the upstream
        pressure p1 are specified and uses an iterative approach with the 
        D-W eqn to solve for the average pressure in the pipe.
        In testing, with typical short pipe segments between orifices it made
        very little difference. For test runs with high flow and long segments
        between orifices it did make a difference. This is the default since
        it doesn't require significantly more computational effort and is more
        broadly applicable.
        '''
        
        self.rho_air = air.rho_air(self.p1,self.air_temp)
        self.Q = air.Q(self.mdot,self.rho_air)
        delta_p = self.pressure_drop
        
        self.p2 = self.p1 - delta_p
        self.p_average = (self.p1 + self.p2)/2
        
        
        #now update with average values
        converged = False
        while converged == False:
            self.rho_air = air.rho_air(self.p_average,self.air_temp)
            self.Q = air.Q(self.mdot,self.rho_air)
            delta_p = self.pressure_drop
            self.p2 = self.p1 - delta_p
            self.p_average_1 = (self.p1 + self.p2)/2
            
            if self.p_average_1 - self.p_average < 0.001:
                 converged = True
            
            self.p_average = self.p_average_1
            
        return delta_p
            
                 
    def __str__(self):
        
        return f'Pipe with ' \
                f'length = {self.length}, ' \
                f'diameter = {self.diameter}, ' \
                f'roughness = {self.roughness}'
    
    @property            
    def f_BUB300_legacy(self):
        f = 0.3164/self.Re**0.25
        k = self.epsilon
        for i in range(50):
            if k == 0:
                k = 0
                f1 = 2* math.log10(self.Re*math.sqrt(f)-0.8)
            else:
                rok = self.diameter/k/2
                f1 = 1.74 + 2*math.log10(rok) - 2*math.log10(1 + 18.7*rok/self.Re/math.sqrt(f))
            if f1 < 10**-20:
                print('f1 < 10e-20')
                return f1
            f1 = 1/f1**2
            if abs(1-f1/f) < 0.05:
                return f1
            f = f1
            