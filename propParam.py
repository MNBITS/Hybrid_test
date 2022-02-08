
from Inputs import UserInputs

grain = UserInputs(10, 10, 10, 60, 30,5100)
grain.setChamberPressure()



class FuelData():
    def __init__(self, Go, c_star, exp, density, coeff):
        self.Go = Go
        self.c_star = c_star
        self.mf_dot = Go*grain.chamberP*grain.At/c_star
        self.coeff = coeff
        self.exp = exp
        self.density = density
        self.thrust_user=thrust_user

    

