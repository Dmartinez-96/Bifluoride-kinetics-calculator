# Biflouride Kinetics Calculator
This calculator operates on a user-provided CSV file to model the kinetics of the thermal decomposition of a bifluoride salt (e.g. NaHF2) into a fluoride salt and hydrogen fluoride gas. There are 15 possible models to choose from in the calculator (see below for info). The source code can be compiled in a Python3 IDE (if preferred), but the calculator itself is a standalone executable that can be run from a bash shell.

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
