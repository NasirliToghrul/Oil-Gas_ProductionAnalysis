import math
import matplotlib.pyplot as plt
import pvtprops as pvt

#Input: Values are the examples values, should be changed according to the requirements

a = 1
b = 2
c = 3
d = 4
e = 5
f = 6
g = 7

tubing_ID = 1.995#inches
wellhead_P = 100 + (g+1)*20 #psi
liquid_rate = 1000 + (d+1)*200 #bbl/d
GLR = 150 #scf/stb
WC = 10 #%
oil_gravity = 20 + (e+1)*2 #API
water_gravity = 1.05
gas_gravity = 0.68
water_gravity = 0
N2=0
CO2=0
H2S=0
water_fvf = 1.02 #rb/stb
tubing_depth = 4000 + (f+1)*500 #ft
Twh = 80 #F
T_gradient = 13/1000

number_of_segments = 5000
segment_length = tubing_depth/number_of_segments
gas_production_rate = GLR*liquid_rate
water_production_rate = (WC/100)*liquid_rate
oil_production_rate = liquid_rate - water_production_rate
producing_GOR =gas_production_rate/oil_production_rate
sg_oil = pvt.calculate_sg_oil(oil_gravity)
WOR = water_production_rate/oil_production_rate
M = 350.17*(sg_oil+WOR*water_gravity) + producing_GOR*gas_gravity*0.0765
numerator_reynolds = 1.4737*10**-5*M*oil_production_rate/(tubing_ID/12)
friction_factor = 4*10**(1.444-2.5*math.log(numerator_reynolds))
K = friction_factor*oil_production_rate**2*M**2/7.4137/10000000000/(tubing_ID/12)**5
wellhead_GOR = pvt.calculate_GOR(producing_GOR,gas_gravity,wellhead_P,oil_gravity,Twh)
wellhead_Bo = pvt.calculate_FVF(producing_GOR,gas_gravity,wellhead_P,oil_gravity, Twh)
wellhead_Z = pvt.redlich_kwong_eos(wellhead_P, Twh, gas_gravity, CO2, H2S, N2)
wellhead_Vm = 5.615*(wellhead_Bo+WOR*water_fvf) + (producing_GOR-wellhead_GOR)*(14.7/wellhead_P)*((Twh+460)/520)*wellhead_Z
wellhead_density = M/wellhead_Vm


pressure_list = [wellhead_P]
temperature_list = [Twh]
GOR_list = [wellhead_GOR]
density_list = [wellhead_density]
depth_list = [0]
average_density_list = []
method = "CarrKobayashi"

depth = 0
previous_density = wellhead_density
P_prime = wellhead_P
T_prime = Twh + segment_length*T_gradient
depth += segment_length

while depth<=tubing_depth:
    while True:
        GOR_prime = pvt.calculate_GOR(producing_GOR,gas_gravity,P_prime,oil_gravity,T_prime)
        compressibility = pvt.redlich_kwong_eos(P_prime, T_prime, gas_gravity, CO2, H2S, N2)
        FVF_oil = pvt.calculate_FVF(producing_GOR, gas_gravity, P_prime, oil_gravity, T_prime)
        V = 5.615*(FVF_oil + WOR*water_fvf) + (producing_GOR - GOR_prime)*(14.7/P_prime)*((T_prime+460)/520)*compressibility
        density = M/V
        density_average = (density+previous_density)/2
        deltaP = (density_average + (K/density_average))*(segment_length/144)
        P_estimated = P_prime + deltaP
        error = abs((P_estimated-P_prime)/P_estimated)
        if error <=0.001:
            pressure_list.append(P_estimated)
            P_prime = P_estimated
            depth_list.append(depth)
            depth += segment_length
            density_list.append(density)
            previous_density = density
            temperature_list.append(T_prime)
            T_prime += segment_length*T_gradient
            GOR_list.append(GOR_prime)
            average_density_list.append(density_average)
            break
        else:
            P_prime = P_estimated


            
plt.plot(pressure_list, depth_list)
plt.title("Pressure Profile")
plt.xlabel("Pressure, psia")
plt.ylabel("Depth, feets")
plt.gca().invert_yaxis()
plt.show()


plt.plot(average_density_list, depth_list[1:])
plt.title("Average Density Distribution")
plt.xlabel("Average Density")
plt.ylabel("Depth")
plt.gca().invert_yaxis()
plt.show()

plt.plot(GOR_list, depth_list)
plt.title("Gas In Solution Profile")
plt.xlabel("Gas In Solution")
plt.ylabel("Depth")
plt.gca().invert_yaxis()
plt.show()

