import numpy as np

class UserInputs():

    def __init__(self, outerDia, coreDia, length, throatarea, stepSize,Thrust_req,total_steps,useBarlow,g,fos,stepNum,Z):
        self.outerDia = outerDia
        self.coreDia = coreDia
        self.length = length
        self.stepSize = stepSize
        self.At = throatarea
        self.portArea = np.pi*(pow(self.outerDia, 2) - pow(self.coreDia, 2))*0.25*self.length
        self.Thrust_req=Thrust_req
        self.total_steps=total_steps
        self.useBarlow=useBarlow
        self.g=g
        self.fos=fos
        self.stepNum=stepNum
        self.Z=Z

    
        
        
        # Grain Volume pi(D2 - d2)/4*length
        self.grainVol = np.pi*((pow(self.outerDia, 2) - pow(self.coreDia, 2)))*0.25*self.length
        


    def setChamberPressure(self):
        self.chamberP = 2.068e+6






