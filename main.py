from Inputs import UserInputs
from propParam import FuelData
import numpy as np
import matplotlib.pyplot as mp
import pandas as pd


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


fuel = FuelData(exp, density, coeff, combustionTemp, k)      # Go, c_star, exponent(burn rate), coefficient, burn rate

#Throat Area Calculation 
throatArea = (Thrust_req)/(fuel.calculate_Cf()*Pc)

print("Throat area: ", throatArea)

grain = UserInputs(outerDia, coreDia, length, throatArea, 0)     # Outer dia, Core Dia, length, throat area, step size


# Grain Volume pi(D2 - d2)/4*length
grainVol = np.pi*((pow(grain.outerDia, 2) - pow(grain.coreDia, 2)))*0.25*grain.length

#Port to throat
Va = np.pi*(grain.outerDia**2)*grain.length*0.25 - grainVol
V1 = grainVol/Va
portTothroat = np.pi*grain.outerDia**2*(1-V1)/(4*throatArea)


#Grain Mass grainVol*Density
fuelMass = grainVol*fuel.density


Isp = 30
OF_init = 8
g = 9.81
c_star = fuel.calculateC_star()

print("c_star: ", c_star)
portArea = np.pi*0.25*pow(grain.coreDia, 2)
print("port area: ", portArea)

m_t = (g*Pc*portArea)/(c_star)
print("mt: ", m_t)
m_ox = (m_t)*OF_init/(OF_init + 1)
Go_init = (m_ox)/(0.25*np.pi*pow(grain.coreDia, 2))
regRate_init = a*pow(Go_init, exp)
print("regrate: ", regRate_init)

BurnTime = (grain.outerDia - grain.coreDia)/(2*regRate_init)

totalSteps = 50
stepSize = BurnTime/totalSteps

print("Burn Time: ", BurnTime)



#initialize Time array
time = np.zeros(totalSteps, dtype = float)

#initiate arrays
radius = np.empty(totalSteps, dtype = float)#radius of port
m_dot_ox = np.zeros(totalSteps, dtype = float)#oxidizer mass flow
m_dot_f = np.zeros(totalSteps, dtype = float)#fuel mass flow
m_dot_t = np.zeros(totalSteps, dtype = float)#total mass flow   
Ph = np.zeros(totalSteps, dtype = float)#
OF_ratio = np.zeros(totalSteps, dtype = float)
regRate = np.zeros(totalSteps, dtype = float)
Gox = np.zeros(totalSteps, dtype = float)
Ptk = np.zeros(totalSteps, dtype = float)
throat_erosion = np.zeros(totalSteps, dtype = float)
thrust = np.zeros(totalSteps,dtype = float)


#OF_desired = float(8)
Z = 1 #Resistence in oxidiser path

#for#OF_ratio = OF_desired

#Calculate Initial Values at t = 0


#" " "  "
m_dot_t[0]=m_t
OF_ratio[0]=OF_init
m_dot_ox[0]=m_ox
m_dot_f[0]=m_t/(1+OF_init)

stepNum = 0
while(stepNum < totalSteps):
    #df = pd.DataFrame({"Thrust" : thrust, "Gox" : Gox,"m_dot_f":m_dot_f,"_dot_ox":m_dot_ox,"m_dot_f":m_dot_f,"m_dot_t":m_dot_t,"OF_ratio":OF_ratio
                       #,"Ph":Ph,"radius":radius,"regRate":regRate,"throat_erosion":throat_erosion,"time":time})
    #df.to_csv("D:\Desktop\Propulsion\Hybrid_sim(csv).csv", index=False)
    
    #Calculate time
    time[stepNum+1] = stepSize*stepNum
    #print("time:" , time)
     
    #Calculate m_dot_t
    m_dot_t[stepNum+1] = Pc*grain.At/(fuel.calculateC_star())
   
 
    #Calculate Ph
    Ph[stepNum] = (1+0.5*(pow(fuel.k*pow(2/fuel.k+1, fuel.k+1/fuel.k-1) , 0.5)*(1/portTothroat)**2))*Pc

    #Calculate M_dotf
    d = 1+OF_ratio[stepNum]
    
    m_dot_f[stepNum+1]=m_dot_t[stepNum+1]/d

    #Calculate Mass flow rate of Oxidiser
    m_dot_ox[stepNum+1] = OF_ratio[stepNum]*m_dot_f[stepNum+1]

    #Calculate G_ox
    Gox[stepNum+1] = m_dot_ox[stepNum]/grain.portArea

    #Calculate regression rate
    regRate[stepNum] = a*pow(Gox[stepNum], fuel.exp)
    
    
    #Calculate burn time
    if (stepNum == 0):
        BurnTime = 0.5*(grain.outerDia - grain.coreDia)/(regRate[stepNum])
    
    #Calculate Pressure of tank
    Ptk = pow(m_dot_ox[stepNum], 2)*Z + Ph[stepNum]

    #Calculate radius
    radius[stepNum] = pow(a*(2*fuel.exp + 1)*pow(m_dot_ox[stepNum]/np.pi , fuel.exp)*time[stepNum] + pow(coreDia/2, 2*fuel.exp + 1), 1/(2*fuel.exp +1))
    throat_erosion[stepNum]=(radius[stepNum]-radius[stepNum-1])
    #Calculate instantaneous mixture ratio
    OF_ratio[stepNum+1] = 1/(2*fuel.density*grain.length*a)*pow(m_dot_ox[stepNum]/np.pi , 1-fuel.exp)*pow(radius[stepNum], 2*fuel.exp - 1)

    #Thrust Calculation
    thrust[stepNum] = m_dot_t[stepNum+1]*Isp*Gox[stepNum]

    #Calculate Port Area For  next iteration
    portArea = np.pi*(pow(grain.outerDia, 2) - pow(4*radius[stepNum], 2))*0.25
    stepNum+=1
    massEjected = m_dot_t[0]*grain.stepSize*stepNum
    fuelMass = np.pi*fuel.density*grain.length*(pow(grain.outerDia, 2)/4 - pow(radius[stepNum], 2))

    print("\n")
  

print(radius)
print(throat_erosion)

# Draw Plots
mp.plot(radius, time, color = "red")
mp.show()

#mp.plot(thrust, time, color = "red")
#mp.show()

mp.plot(OF_ratio, time, color = "red")
mp.show()

mp.plot(regRate, time, color = "red")

mp.show()
mp.show()
