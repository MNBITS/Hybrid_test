
from Inputs import UserInputs

grain = UserInputs(10, 10, 10, 60, 30)
grain.setChamberPressure()
R = 8.31


class FuelData():

    def __init__(self, Go, exp, density, coeff, combustionTemp,k):

        self.Go = Go
        #self.mf_dot = Go*grain.chamberP*grain.At/c_star
        self.coeff = coeff
        self.exp = exp
        self.density = density
        self.combustionTemp = combustionTemp
        self.k = k


    def calculateC_star(self):
        numerator = self.combustionTemp*R
        alpha = self.k+1/self.k-1
        beta = 2/self.k+1


        c_star = pow(numerator/self.k*pow(beta, alpha), 0.5)

        return c_star

    def calculate_Cf(self):
        alpha = self.k+1/self.k-1
        print("alpha: ", alpha)
        beta = 2/self.k+1
        gamma = 2*pow(self.k, 2)/self.k - 1
        eta = pow((100000)/(grain.chamberP), (self.k - 1)/(self.k))
        print("eta: ", eta)
        #Cf=((2*self.k**2/(self.k-1))*((beta)**(alpha))*(1-(100000/grain.chamberP)**((self.k-1)/self.k)))**0.5

        Cf = pow(gamma*pow(beta, alpha)*(1-eta), 0.5)
        print("Cf: ", Cf)
    
        return Cf
       
    

