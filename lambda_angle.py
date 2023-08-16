import numpy as np
import matplotlib.pyplot as plt
from utils.analysis import Measurements, Lambda_max, Angles

'''
code to find optimal value of the prism's angle from experimental data
as a wavelength function

we study several maximums for a given external angle (24° in this case)
and try to minimice the difference between all the distances.

inestable method to find the optimal angle.
'''

# read experimental data 

lambda_T_dict = Measurements('data/sources')._lambda_transmittance()

wavelength, T_exp = lambda_T_dict['24.0']

plt.figure()
plt.plot(wavelength, T_exp, '.')
plt.grid(alpha = 0.7)
plt.title(r'Experimental transmittance for $\theta_{ext}$ = 24.0°' )
plt.xlabel(r'$\lambda (nm)$')
plt.xlim(440,1000)
plt.ylabel('T')
plt.show()


# we obtain the values of the wavelength for each maximum (5 in total)

values = [450, 500, 560, 650, 780, 920] #values of lambda that delimiter the maximums
maximums =[]

for i in range(len(values)):
    try:
        _ , _ , lambda_1 , lambda_2 = Lambda_max( T = T_exp, wavelength= wavelength).parabolic(x = wavelength,
                                                                                              value_1 = values[i],
                                                                                              value_2 = values[i+1],
                                                                                              value_3 = values[i+2])

        print(lambda_1, lambda_2)

        if round(lambda_1*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_1*1e-3)

        if round(lambda_2*1e-3, 3) not in [round(val, 3) for val in maximums]: 
            maximums.append(lambda_2*1e-3)
    
    except Exception:
        break

maximums = maximums[1:-1] # we don't consider first and last maximums for being too close to the limits


'''
now we compute for each value of alpha the distances and internal angles. 

the optimal value of alpha will be the one that minimices the
standar deviation for the d's.
'''


a = np.linspace(44.2,50,100)*np.pi/180
theta_ext = 24 * np.pi/180
min_0 = 1000

for alpha in a:
    print('\n', alpha*180/np.pi)
    d = []

    for i in range(len(maximums)):
        try:
            theta_1 = Angles(alpha, wavelength = maximums[i], theta_ext = theta_ext).int_angle()
            theta_2 = Angles(alpha, wavelength = maximums[i+1], theta_ext = theta_ext).int_angle()

            dist = Lambda_max( T = T_exp, wavelength = wavelength).distance(theta_1=theta_1, theta_2=theta_2,
                                                                           lambda_1 = maximums[i], lambda_2 = maximums[i+1])
            
            d.append(dist)
            # print(d)

        except Exception:
            pass

    d_array = np.array(d)
    d_mean = np.mean(d_array)

    min = np.sqrt(np.sum((d_array-d_mean)**2)/len(d))
    # print('min', min)
    # print('min0', min_0)

    if min < min_0:
        final_a = alpha
        final_d = d
        min_0 = min
    else:
        pass

print(final_a*180/np.pi)
print(final_d)
