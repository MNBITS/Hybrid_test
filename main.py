from Inputs import UserInputs
from propParam import FuelData
import numpy as np
import matplotlib.pyplot as mp
from random import random

coreDia = 0
outerDia = 10
length = 0
stepSize = 0.
Thrust = 1000
Cf = 1
useBarlow = True
fos = 2
maxStress = 0
Pc = 0

def setChamberPressure(useBarlow, fos):
    if(useBarlow):
        Pc = 2*maxStress*(outerDia - coreDia)/(outerDia*fos)
    
    else:
        Pc = 2000
    return Pc

throatArea = Thrust/Cf*setChamberPressure(20,30)

exp = 0.2
Go = 120
a = 0
burnRate = 0
density = 0
coeff = 0
combustionTemp = 0
k = 0


grain = UserInputs(outerDia, coreDia, length, throatArea, stepSize)     # Outer dia, Core Dia, length, throat area, step size
fuel = FuelData(Go, exp, density, coeff, combustionTemp, k)      # Go, c_star, exponent(burn rate), coefficient, burn rate

# Grain Volume pi(D2 - d2)/4*length
grainVol = np.pi*((pow(grain.outerDia, 2) - pow(grain.coreDia, 2)))*0.25*grain.length

#Port to throat
Va = np.pi*grain.outerDia*grain.length - grainVol
V1 = grainVol/Va
portTothroat = np.pi*grain.outerDia**2*(1-V1)/(4*throatArea)


#Port to throat
Va = np.pi*grain.outerDia*grain.length - grainVol
V1 = grainVol/Va
portTothroat = np.pi*grain.outerDia**2*(1-V1)/(4*throatArea)


# Grain Mass grainVol*Density
fuelMass = grainVol*fuel.density

m_dot = 90
Isp = 0
BurnTime = 20

totalSteps = int(BurnTime/stepSize)

#initialize Time array
time = np.zeros(totalSteps, dtype = float)

stepNum = 0
while(stepNum < totalSteps):
    
    #initiate arrays
    radius = np.empty(totalSteps, dtype = float)
    m_dot_ox = np.zeros(totalSteps, dtype = float)
    m_dot_f = np.zeros(totalSteps, dtype = float)
    m_dot_t = np.zeros(totalSteps, dtype = float)
    Ph = np.zeros(totalSteps, dtype = float)
    OF_ratio = np.zeros(grain.stepSize, dtype = float)
     
    #Calculate m_dot_t
    m_dot_t[stepNum] = fuel.Go*Pc*throatArea/(fuel.calculateC_star)
    
    #Calculate time
    time[stepNum] = grain.stepSize*stepNum

    #Calculate Ph
    Ph[stepNum] = [1+0.5*pow(fuel.k*pow(2/fuel.k+1), fuel.k+1/fuel.k-1)]

    #m_dot_f[stepNum] = 2*np.pi*fuel.density*grain.length*alpha*pow(radius[stepNum], 1-2*fuel.exp)





    #m_dot_tot[stepNum] = m_dot_f[stepNum] + m_dot_ox[stepNum]


    thrust = np.zeros(grain.stepSize, dtype = float)

    thrust[stepNum] = m_dot_t[stepNum]*Isp*fuel.Go



    portArea = np.pi*(pow(grain.outerDia, 2) - pow(4*radius[stepNum], 2))*0.25
    stepNum+=1
    massEjected = m_dot*grain.stepSize*stepNum
    fuelMass = np.pi*fuel.density*grain.length*(pow(radius[stepNum], 2) - pow(grain.coreDia, 2))
    burnRate=m_dot_t[stepNum]/((throatArea)*density)
    
    #print(radius[stepNum])
    print("\n")

mp.plot(thrust, time, color = "red")
mp.show()
