import math
import matplotlib.pyplot as plt
import pvtprops as pvt  # Importing custom module for PVT properties calculations

# Input: Values are example values, should be changed according to the requirements

tubing_ID = 1.995  # inches
wellhead_P = 260  # psi
liquid_rate = 2000  # bbl/d
GLR = 150  # scf/stb (Gas-Liquid Ratio)
WC = 10  # % Water Cut
oil_gravity = 32  # API (American Petroleum Institute)
water_gravity = 1.05  # Water Specific Gravity
gas_gravity = 0.68  # Gas Specific Gravity
N2 = 0  # Nitrogen mole fraction (for PVT calculations)
CO2 = 0  # Carbon Dioxide mole fraction (for PVT calculations)
H2S = 0  # Hydrogen Sulfide mole fraction (for PVT calculations)
water_fvf = 1.02  # rb/stb (Reservoir barrels per stock tank barrel)
tubing_depth = 7500  # ft
Twh = 80  # F (Wellhead temperature)
T_gradient = 13 / 1000  # Temperature gradient per foot

number_of_segments = 5000
segment_length = tubing_depth / number_of_segments  # ft
gas_production_rate = GLR * liquid_rate  # scf/d (Standard cubic feet per day)
water_production_rate = (WC / 100) * liquid_rate  # bbl/d
oil_production_rate = liquid_rate - water_production_rate  # bbl/d
producing_GOR = gas_production_rate / oil_production_rate  # Gas-Oil Ratio
sg_oil = pvt.calculate_sg_oil(oil_gravity)  # Specific gravity of oil
WOR = water_production_rate / oil_production_rate  # Water-Oil Ratio
M = 350.17 * (sg_oil + WOR * water_gravity) + producing_GOR * gas_gravity * 0.0765  # Mixture density constant (lbm/ft^3)
numerator_reynolds = 1.4737 * 10 ** -5 * M * oil_production_rate / (tubing_ID / 12)  # Reynolds number numerator
friction_factor = 4 * 10 ** (1.444 - 2.5 * math.log(numerator_reynolds))  # Friction factor
K = friction_factor * oil_production_rate ** 2 * M ** 2 / 7.4137 / 10000000000 / (tubing_ID / 12) ** 5  # K parameter
wellhead_GOR = pvt.calculate_GOR(producing_GOR, gas_gravity, wellhead_P, oil_gravity, Twh)  # Wellhead Gas-Oil Ratio
wellhead_Bo = pvt.calculate_FVF(producing_GOR, gas_gravity, wellhead_P, oil_gravity, Twh)  # Wellhead Formation Volume Factor
wellhead_Z = pvt.redlich_kwong_eos(wellhead_P, Twh, gas_gravity, CO2, H2S, N2)  # Wellhead Z factor
wellhead_Vm = 5.615 * (wellhead_Bo + WOR * water_fvf) + (producing_GOR - wellhead_GOR) * (14.7 / wellhead_P) * (
            (Twh + 460) / 520) * wellhead_Z  # Wellhead mixture volume
wellhead_density = M / wellhead_Vm  # Wellhead fluid density (lbm/ft^3)

# Initializing lists to store pressure, temperature, GOR, density, and depth values
pressure_list = [wellhead_P]  # psi
temperature_list = [Twh]  # F
GOR_list = [wellhead_GOR]  # Gas-Oil Ratio
density_list = [wellhead_density]  # lbm/ft^3
depth_list = [0]  # ft
average_density_list = []  # lbm/ft^3
method = "CarrKobayashi"  # Calculation method (unused in this script)

depth = 0
previous_density = wellhead_density
P_prime = wellhead_P
T_prime = Twh + segment_length * T_gradient
depth += segment_length

# Iterating through depth segments to calculate pressure profiles
while depth <= tubing_depth:
    while True:
        # Calculating parameters at current depth
        GOR_prime = pvt.calculate_GOR(producing_GOR, gas_gravity, P_prime, oil_gravity, T_prime)
        compressibility = pvt.redlich_kwong_eos(P_prime, T_prime, gas_gravity, CO2, H2S, N2)
        FVF_oil = pvt.calculate_FVF(producing_GOR, gas_gravity, P_prime, oil_gravity, T_prime)
        V = 5.615 * (FVF_oil + WOR * water_fvf) + (producing_GOR - GOR_prime) * (14.7 / P_prime) * (
            (T_prime + 460) / 520) * compressibility  # Mixture volume
        density = M / V  # Fluid density
        density_average = (density + previous_density) / 2  # Average density
        deltaP = (density_average + (K / density_average)) * (segment_length / 144)  # Pressure change
        P_estimated = P_prime + deltaP  # Estimated pressure
        error = abs((P_estimated - P_prime) / P_estimated)  # Error calculation
        if error <= 0.001:
            # Storing calculated values
            pressure_list.append(P_estimated)
            P_prime = P_estimated
            depth_list.append(depth)
            depth += segment_length
            density_list.append(density)
            previous_density = density
            temperature_list.append(T_prime)
            T_prime += segment_length * T_gradient
            GOR_list.append(GOR_prime)
            average_density_list.append(density_average)
            break
        else:
            P_prime = P_estimated

# Plotting pressure profile
plt.plot(pressure_list, depth_list)
plt.title("Pressure Profile")
plt.xlabel("Pressure (psia)")
plt.ylabel("Depth (feet)")
plt.gca().invert_yaxis()
plt.show()

# Plotting average density distribution
plt.plot(average_density_list, depth_list[1:])
plt.title("Average Density Distribution")
plt.xlabel("Average Density (lbm/ft^3)")
plt.ylabel("Depth (feet)")
plt.gca().invert_yaxis()
plt.show()

# Plotting Gas In Solution profile
plt.plot(GOR_list, depth_list)
plt.title("Gas In Solution Profile")
plt.xlabel("Gas In Solution (scf/stb)")
plt.ylabel("Depth (feet)")
plt.gca().invert_yaxis()
plt.show()
