# Bifluoride Kinetics Calculator
This calculator operates on a user-provided CSV file to model the kinetics of the thermal decomposition of a bifluoride salt (e.g. NaHF2) into a fluoride salt and hydrogen fluoride gas as a best fit curve. There are 15 possible models to choose from in the calculator (see below for info) for the sake of comparison. The source code can be compiled in a Python3 IDE (if preferred), but the calculator itself is a standalone executable that can be run from a bash shell.

## Calculator Scope:
The calculator is designed to read in a CSV file of a specific format containing experimental data on the thermal decomposition of a given bifluoride salt. 

The first four columns of your CSV file need to be labeled exactly as follows (see example below):

| time_data | initmol_data | temp_data | finalmol_data |
| ------ | ------ | ------ | ------ |

Your CSV file should contain your experimental data in the order mentioned above. The respective units your data should be in are: 

| Time taken to run full experiment in seconds | Initial moles of bifluoride salt used | Constant experimental temperature in Kelvin | Final moles of HF evolved in experiment |
| ------ | ------ | ------ | ------ |

Your CSV file should look similar to the below image (with your own data input):

![image](https://user-images.githubusercontent.com/85904612/123525859-808d0d00-d699-11eb-8c80-85f96124fae3.png)

The program will output the best-fit parameters for the reaction activation energy (in kcal/mol), Ea, and the Arrhenius pre-factor (dimensionless), A. Additionally, the output includes the variance of the parameters _A_ and _Ea_, as well as the covariance between _A_ and _Ea_. These metrics can be used to determine how good of a fit you have obtained with your input data and selected model.

The calculator models reaction progress, parameterized by alpha, as a function of how long the experiment was run (time in seconds) and constant experimental temperature (in Kelvin), based on the different models described in the section "Models Available" below.

The produced plot models your data from your data CSV file, scaled such that all data points are at the average temperature of the data points in your file. 

## Package Installation for Running bifluoride_kinetics_curve_fitting_source_code.py in Python3 IDE
The source code can be compiled in a Python 3 IDE. The required packages are:

```sh
numpy
scipy
warnings
pandas
matplotlib
csv
```

The easiest way to install these packages is to use the pip install comand:

```sh
pip install <package>
```

## Instructions for Running bifluoride_kinetics_curve_fitting_source_code.py in Python3 IDE
- Obtain a CSV file of your experimental data in the format mentioned above, in the Calculator Scope section. 
- Compile the source code in your IDE. You will be prompted to enter the directory containing your data CSV file. See example below:

![image](https://user-images.githubusercontent.com/85904612/123526538-30fd1000-d69e-11eb-805f-16f6f00a7d4d.png)

- Use the full directory if your CSV file is outside of your Python path. See example below:

![image](https://user-images.githubusercontent.com/85904612/123526592-8e915c80-d69e-11eb-8634-c4572c9c0071.png)

- The calculator will then prompt you to select the model number you would like to test with your data. See example below:

![image](https://user-images.githubusercontent.com/85904612/123526632-d4e6bb80-d69e-11eb-9bae-8d06be7f7863.png)

- Choose a particular model and press enter. If you entered a valid number, the calculator will produce your best fit parameter values as well as informing you that the resultant plot for the selected model will be saved in the current working directory. You will also be prompted to try a new model with the same data, use a new data file, or quit. See example below:

![image](https://user-images.githubusercontent.com/85904612/123526722-65250080-d69f-11eb-97af-b4eeb45cd63a.png)

Example plot produced: 

![Zero-order_CurveFit](https://user-images.githubusercontent.com/85904612/123526783-b9c87b80-d69f-11eb-8b6d-bb842d2bb9fa.png)

- If you input 'Y' or 'y' without the quotations, you will be returned to the model selection screen with the same data file. If you input 'C' or 'c' without the quotations, you will again be prompted to input the directory of your data CSV file. If you input 'N' or 'n' the program will quit. If you input anything else, the program will simply return you to the model selection screen.

## Instructions for Running bifluoride_kinetics_curve_fitting in Terminal:
- Module dependencies are already packaged into bifluoride_kinetics_curve_fitting.
- Obtain a CSV file of your experimental data in the format mentioned above, in the Calculator Scope section. 
- Run bash or open up terminal. See example below:

![image](https://user-images.githubusercontent.com/85904612/123526933-dc0ec900-d6a0-11eb-8e8b-54f61395afa0.png)

- cd to location of bifluoride_kinetics_curve_fitting file. See example below:

![image](https://user-images.githubusercontent.com/85904612/123526971-2728dc00-d6a1-11eb-8b4d-3920e1b223db.png)

- Run the following command:

```sh
./bifluoride_kinetics_curve_fitting
```

- Input the directory of your data CSV file when prompted. See example below:

![image](https://user-images.githubusercontent.com/85904612/123527021-99012580-d6a1-11eb-8b96-826e7a39cf35.png)

- The calculator will then prompt you to select the model number you would like to test with your data. See example below:

![image](https://user-images.githubusercontent.com/85904612/123527073-217fc600-d6a2-11eb-8d83-dc124d1bfe24.png)

- Choose a particular model and press enter. If you entered a valid number, the calculator will produce your best fit parameter values as well as informing you that the resultant plot for the selected model will be saved in the current working directory. You will also be prompted to try a new model with the same data, use a new data file, or quit. See example below:

![image](https://user-images.githubusercontent.com/85904612/123527108-6b68ac00-d6a2-11eb-8bf4-e9e01ecf566f.png)

- If you input 'Y' or 'y' without the quotations, you will be returned to the model selection screen with the same data file. If you input 'C' or 'c' without the quotations, you will again be prompted to input the directory of your data CSV file. If you input 'N' or 'n' the program will quit. If you input anything else, the program will simply return you to the model selection screen.

## Models Available:
There are 15 models to select from in this program. They are each identified by a particular ID number, mentioned in the instructions above. Below is a table specifying which ID number corresponds to which model:

| ID Number | Kinetics Model Name | Model Function |
| ------ | ------ | ------ |
| 1 | Zero-order | ![image](https://user-images.githubusercontent.com/85904612/123526006-a0710080-d69a-11eb-836f-7322570d5072.png) |
| 2 | First-order | ![image](https://user-images.githubusercontent.com/85904612/123526051-f9d92f80-d69a-11eb-95dd-6c2d01482348.png) |
| 3 | Second-order | ![image](https://user-images.githubusercontent.com/85904612/123526057-02ca0100-d69b-11eb-9036-70040f4bf1df.png) |
| 4 | Third-order | ![image](https://user-images.githubusercontent.com/85904612/123526086-28efa100-d69b-11eb-931a-54a050867b30.png) |
| 5 | Avrami-Erofeyev 1 | ![image](https://user-images.githubusercontent.com/85904612/123526096-32790900-d69b-11eb-8af6-a8c9d9231bc4.png) |
| 6 | Avrami-Erofeyev 2 | ![image](https://user-images.githubusercontent.com/85904612/123526101-373dbd00-d69b-11eb-834e-9520635cc0d7.png) |
| 7 | Avrami-Erofeyev 3 | ![image](https://user-images.githubusercontent.com/85904612/123526180-b3d09b80-d69b-11eb-9492-ea57890564c8.png) |
| 8 | Avrami-Erofeyev 4 | ![image](https://user-images.githubusercontent.com/85904612/123526118-59373f80-d69b-11eb-8b25-c03669333e11.png) |
| 9 | 2/3 power law | ![image](https://user-images.githubusercontent.com/85904612/123526138-74a24a80-d69b-11eb-9616-a356b7618eb4.png) |
| 10 | Quadratic power law | ![image](https://user-images.githubusercontent.com/85904612/123526143-8257d000-d69b-11eb-8dc9-ab7394fca463.png) |
| 11 | Cubic power law | ![image](https://user-images.githubusercontent.com/85904612/123526158-93a0dc80-d69b-11eb-80a0-e1bcb269e523.png) |
| 12 | Quartic power law | ![image](https://user-images.githubusercontent.com/85904612/123526167-9dc2db00-d69b-11eb-9f83-351debf3231a.png) |
| 13 | Contracting area | ![image](https://user-images.githubusercontent.com/85904612/123526189-c34fe480-d69b-11eb-8b6a-c848add3430e.png) |
| 14 | Contracting volume | ![image](https://user-images.githubusercontent.com/85904612/123526194-ca76f280-d69b-11eb-82b6-ddd00e708b16.png) |
| 15 | 1-D diffusion | ![image](https://user-images.githubusercontent.com/85904612/123526199-d1056a00-d69b-11eb-83e7-47d9c472a08b.png) |
