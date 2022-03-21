from Inputs import UserInputs
from propParam import FuelData
import numpy as np
import matplotlib.pyplot as mp


coreDia = 0.052
outerDia = 0.067
length = 0.6



Thrust_req = 1000
useBarlow = True
fos = 2
maxStress = 230
Pc = 2.068000e+006


def setChamberPressure(useBarlow, fos):
    if(useBarlow):
        Pc = 2*maxStress*(outerDia - coreDia)/(outerDia*fos)
    
    else:
        Pc = 2.068e+6
    return Pc



exp = 0.62
#Go = 120
a = 0.117
burnRate = 30
density = 788.6
coeff = 20
combustionTemp = 300
k = 1.250
Isp = 30
OF_init = 8
stepNum = 0
Z = 1 #Resistence in oxidiser path
g = 9.81
totalSteps = 50


fuel = FuelData(exp, density, coeff, combustionTemp, k,burnRate,a,Isp,OF_init)      # Go, c_star, exponent(burn rate), coefficient, burn rate

#Throat Area Calculation 
throatArea = (Thrust_req)/(fuel.calculate_Cf()*Pc)

print("Throat area: ", throatArea)

grain = UserInputs(outerDia, coreDia, length, throatArea, 0,Thrust_req,totalSteps,useBarlow,g,fos,stepNum,Z)     # Outer dia, Core Dia, length, throat area, step size


# Grain Volume pi(D2 - d2)/4*length
grainVol = np.pi*((pow(grain.outerDia, 2) - pow(grain.coreDia, 2)))*0.25*grain.length

#Port to throat
Va = np.pi*(grain.outerDia**2)*grain.length*0.25 - grainVol
V1 = grainVol/Va
portTothroat = np.pi*grain.outerDia**2*(1-V1)/(4*throatArea)


#Grain Mass grainVol*Density
fuelMass = grainVol*fuel.density





c_star = fuel.calculateC_star()

print("c_star: ", c_star)
portArea = np.pi*0.25*pow(grain.coreDia, 2)
print("port area: ", portArea)

m_t = (g*Pc*throatArea)/(c_star)
print("mt: ", m_t)
m_ox = (m_t)/(OF_init + 1)
Go_init = (m_ox)/(0.25*np.pi*pow(grain.coreDia, 2))
regRate_init = a*pow(Go_init, exp)
print("regrate: ", regRate_init)

BurnTime = (grain.outerDia - grain.coreDia)/(2*regRate_init)


stepSize = BurnTime/totalSteps

print("Burn Time: ", BurnTime)



#initialize Time array
time = np.zeros(totalSteps, dtype = float)

#initiate arrays
radius = np.empty(totalSteps, dtype = float)
m_dot_ox = np.zeros(totalSteps, dtype = float)
m_dot_f = np.zeros(totalSteps, dtype = float)
m_dot_t = np.zeros(totalSteps, dtype = float)
Ph = np.zeros(totalSteps, dtype = float)
OF_ratio = np.zeros(totalSteps, dtype = float)
regRate = np.zeros(totalSteps, dtype = float)
Gox = np.zeros(totalSteps, dtype = float)
Ptk = np.zeros(totalSteps, dtype = float)
thrust = np.zeros(grain.stepSize, dtype = float)


#OF_desired = float(8)


#for#OF_ratio = OF_desired

#Calculate Initial Values at t = 0

#" " "  "


while(stepNum < totalSteps):
    
    #Calculate time
    time[stepNum] = stepSize*stepNum
    #print("time:" , time)
     
    #Calculate m_dot_t
    m_dot_t[stepNum] = Pc*grain.At/(fuel.calculateC_star())

    #Calculate Ph
    Ph[stepNum] = (1+0.5*(pow(fuel.k*pow(2/fuel.k+1, fuel.k+1/fuel.k-1) , 0.5)*(1/portTothroat)**2))*Pc

    #Calculate M_dotf
    d = 1+OF_ratio[stepNum]


    #Calculate Mass flow rate of Oxidiser
    m_dot_ox[stepNum] = OF_ratio[stepNum]*m_dot_f[stepNum]

    #Calculate G_ox
    Gox[stepNum] = m_dot_ox[stepNum]/grain.portArea

    #Calculate regression rate
    regRate[stepNum] = a*pow(Gox[stepNum], fuel.exp)
    

    #Calculate burn time
    #if (stepNum == 0):
     #   BurnTime = 0.5*(grain.outerDia - grain.coreDia)/(regRate[stepNum])
    
    #Calculate Pressure of tank
    Ptk = pow(m_dot_ox[stepNum], 2)*Z + Ph[stepNum]

    #Calculate radius
    radius[stepNum] = pow(a*(2*fuel.exp + 1)*pow(m_dot_ox[stepNum]/np.pi , fuel.exp)*time[stepNum] + pow(coreDia/2, 2*fuel.exp + 1), 1/(2*fuel.exp +1))

    #Calculate instantaneous mixture ratio
    OF_ratio[stepNum] = 1/(2*fuel.density*grain.length*a)*pow(m_dot_ox[stepNum]/np.pi , 1-fuel.exp)*pow(radius[stepNum], 2*fuel.exp - 1)

    #Thrust Calculation
    #thrust[stepNum] = m_dot_t[stepNum]*Isp*fuel.Go

    #Calculate Port Area For  next iteration
    portArea = np.pi*(pow(grain.outerDia, 2) - pow(4*radius[stepNum], 2))*0.25
    stepNum+=1
    #massEjected = m_dot*grain.stepSize*stepNum
    #fuelMass = np.pi*fuel.density*grain.length*(pow(radius[stepNum], 2) - pow(grain.coreDia, 2))

    #print("\n")
  

print(radius)

# Draw Plots
mp.plot(radius, time, color = "red")
mp.show()

mp.plot(thrust, time, color = "red")
mp.show()

mp.plot(OF_ratio, time, color = "red")
mp.show()

mp.plot(regRate, time, color = "red")

mp.show() 


