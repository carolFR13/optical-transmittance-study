import numpy as np
import matplotlib.pyplot as plt
from optics import Measurements, Lambda_max, Angles, weighted_average
from uncertainties import ufloat


'''
we compute the distances for the lambda mesurements with 
the angle obtained from the He-Ne data.

we assume an incertainty in the external angle of 0.1°
according to the manufacturer.
'''

# read experimental data 

lambda_T_dict = Measurements('data/sources')._lambda_transmittance()

wavelength, T_exp_24 = lambda_T_dict['24.0']
wavelength, T_exp_11 = lambda_T_dict['11.0']
wavelength, T_exp_8 = lambda_T_dict['8.0']


# 24°
# we obtain the values of the wavelength for each maximum (5 in total)

values = [450, 500, 560, 650, 780, 920] #values of lambda that delimiter the maximums
maximums =[]

for i in range(len(values)):
    try:
        _ , _ , lambda_1 , lambda_2 = Lambda_max( T = T_exp_24, wavelength= wavelength).parabolic(x = wavelength,
                                                                                              value_1 = values[i],
                                                                                              value_2 = values[i+1],
                                                                                              value_3 = values[i+2])

        if round(lambda_1*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_1*1e-3)

        if round(lambda_2*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_2*1e-3)
    
    except Exception:
        break

print('24°:',maximums, len(maximums))

alpha = 45.214414414414414 * np.pi/180
theta_ext = ufloat(24, 0.1/(12**(1/2))) * np.pi/180

d = []

for i in range(len(maximums)):
    try:
        theta_1 = Angles(alpha, wavelength = maximums[i], theta_ext = theta_ext).int_angle(uncertainty=True)
        theta_2 = Angles(alpha, wavelength = maximums[i+1], theta_ext = theta_ext).int_angle(uncertainty=True)

        dist = Lambda_max( T = T_exp_24, wavelength = wavelength).distance(theta_1=theta_1, theta_2=theta_2,
                                                                           lambda_1 = maximums[i], lambda_2 = maximums[i+1], uncertainty=True)
            
        d.append(dist)
        # print(d)

    except Exception:
        pass

print('distances:',d)

values = [x.nominal_value for x in d]
std = [x.std_dev for x in d]
mean, u = weighted_average(values, std)

print('mean:', mean, u)


# 11°

values = [450, 520, 700, 1000] #values of lambda that delimiter the maximums
maximums =[]

for i in range(len(values)):
    try:
        _ , _ , lambda_1 , lambda_2 = Lambda_max( T = T_exp_11, wavelength= wavelength).parabolic(x = wavelength,
                                                                                              value_1 = values[i],
                                                                                              value_2 = values[i+1],
                                                                                              value_3 = values[i+2])

        if round(lambda_1*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_1*1e-3)

        if round(lambda_2*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_2*1e-3)
    
    except Exception:
        break
print('11°:',maximums, len(maximums))


theta_ext = ufloat(11, 0.1/(12**(1/2))) * np.pi/180
d = []
for i in range(len(maximums)):
    try:
        theta_1 = Angles(alpha, wavelength = maximums[i], theta_ext = theta_ext).int_angle(uncertainty=True)
        theta_2 = Angles(alpha, wavelength = maximums[i+1], theta_ext = theta_ext).int_angle(uncertainty=True)

        dist = Lambda_max( T = T_exp_11, wavelength = wavelength).distance(theta_1 = theta_1, theta_2 = theta_2,
                                                                           lambda_1 = maximums[i], lambda_2 = maximums[i+1],
                                                                           uncertainty=True)
            
        d.append(dist)
        # print(d)

    except Exception:
        pass

print('distances:',d)

values = [x.nominal_value for x in d]
std = [x.std_dev for x in d]
mean, u = weighted_average(values, std)

print('mean:', mean, u)


# 8°
values = [450, 650, 1000] #values of lambda that delimiter the maximums
maximums =[]

for i in range(len(values)):
    try:
        _ , _ , lambda_1 , lambda_2 = Lambda_max( T = T_exp_8, wavelength= wavelength).parabolic(x = wavelength,
                                                                                              value_1 = values[i],
                                                                                              value_2 = values[i+1],
                                                                                              value_3 = values[i+2])

        if round(lambda_1*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_1*1e-3)

        if round(lambda_2*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_2*1e-3)
    
    except Exception:
        break

print('8°:',maximums, len(maximums))

theta_ext = ufloat(8, 0.1/(12**(1/2))) * np.pi/180
d = []
for i in range(len(maximums)):
    try:
        theta_1 = Angles(alpha, wavelength = maximums[i], theta_ext = theta_ext).int_angle(uncertainty=True)
        theta_2 = Angles(alpha, wavelength = maximums[i+1], theta_ext = theta_ext).int_angle(uncertainty=True)

        dist = Lambda_max( T = T_exp_8, wavelength = wavelength).distance(theta_1=theta_1, theta_2=theta_2,
                                                                           lambda_1 = maximums[i], lambda_2 = maximums[i+1],
                                                                           uncertainty=True)
            
        d.append(dist)
        # print(d)

    except Exception:
        pass

print('ditances:',d)

values = [x.nominal_value for x in d]
std = [x.std_dev for x in d]
mean, u = weighted_average(values, std)

print('mean:', mean, u)