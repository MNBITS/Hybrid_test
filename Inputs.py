import numpy as np

class UserInputs():

    def __init__(self, outerDia, coreDia, length, throatarea,  stepSize):
        self.outerDia = outerDia
        self.coreDia = coreDia
        self.length = length
        self.stepSize = stepSize
        self.At = throatarea
        self.portArea = np.pi*(pow(self.outerDia, 2) - pow(self.coreDia, 2))*0.25*self.length
        


    def setChamberPressure(self):
        self.chamberP = 50





