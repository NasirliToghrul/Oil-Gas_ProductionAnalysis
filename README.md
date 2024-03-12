# Production PVTProps Analysis

This Python script is designed to perform Pressure-Volume-Temperature (PVT) analysis for oil and gas wells using various properties and equations. It calculates pressure, density, and gas in solution profiles along the wellbore depth.

## Table of Contents

- [Requirements](#requirements)
- [Usage](#usage)
- [Input Parameters](#input-parameters)
- [Output](#output)
- [Equations and Methods](#equations-and-methods)

## Requirements

- Python 3.x
- matplotlib library
- pvtprops: A custom library for PVT (Pressure-Volume-Temperature) properties. Details can be found [here](https://github.com/NasirliToghrul/PVTProps_Library).

## Usage

1. Make sure you have Python installed on your system. You can download Python from [python.org](https://www.python.org/downloads/).

2. Install the required libraries using pip:
    ```
    pip install matplotlib
    ```

3. Download or clone the `pvtprops.py` module and ensure it is in the same directory as your script.

4. Run the script by executing the Python file in your terminal or IDE.

    ```
    python production_analysis.py
    ```

## Input Parameters

- **Tubing Inner Diameter (`tubing_ID`)**: Diameter of the tubing in inches.
- **Wellhead Pressure (`wellhead_P`)**: Initial pressure at the wellhead in psi.
- **Liquid Rate (`liquid_rate`)**: Rate of liquid production in barrels per day (bbl/d).
- **Gas-Liquid Ratio (`GLR`)**: Gas-liquid ratio in standard cubic feet per stock tank barrel (scf/stb).
- **Water Cut (`WC`)**: Percentage of water in the produced fluid.
- **Oil Gravity (`oil_gravity`)**: Specific gravity of oil in API gravity.
- **Water Gravity (`water_gravity`)**: Specific gravity of water.
- **Gas Gravity (`gas_gravity`)**: Specific gravity of gas.
- **Water Formation Volume Factor (`water_fvf`)**: Formation volume factor of water in rb/stb.
- **Tubing Depth (`tubing_depth`)**: Total depth of the tubing in feet.
- **Wellhead Temperature (`Twh`)**: Initial temperature at the wellhead in Fahrenheit.
- **Temperature Gradient (`T_gradient`)**: Rate of change of temperature with depth in Fahrenheit per foot.

## Output

The script generates three plots:

1. **Pressure Profile**: Graph showing the variation of pressure with depth along the wellbore.
2. **Average Density Distribution**: Graph illustrating the distribution of average density with depth.
3. **Gas In Solution Profile**: Graph displaying the variation of gas in solution with depth.

## Equations and Methods

The script utilizes various equations and methods to perform the PVT analysis:

- **Redlich-Kwong Equation of State**: Used to calculate the compressibility factor of the gas.
- **Carr-Kobayashi Correlation**: Utilized for estimating pressure drop due to friction in the tubing.
- **Reynolds Number Calculation**: Calculates the Reynolds number for flow regime determination.
- **Formation Volume Factor (FVF)**: Calculates the FVF of oil.
- **Gas-Oil Ratio (GOR) Calculation**: Determines the gas-oil ratio at different depths.
- **Density Calculation**: Calculates fluid density using the ideal gas law.

For more details on the underlying theory and equations, please refer to the comments in the Python script.

#### Additional Notes
- The PVTProps library is utilized for handling fluid properties. Detailed information about the library can be found on [GitHub](https://github.com/NasirliToghrul/PVTProps_Library).

## Example Graphs

![AverageDensityDistribution.png](ExampleGraphs%2FAverageDensityDistribution.png)

![GasInSolutionProfile.png](ExampleGraphs%2FGasInSolutionProfile.png)

![PressureProfile.png](ExampleGraphs%2FPressureProfile.png)

## Credits

This script was developed by Toghrul Nasirli.

---

