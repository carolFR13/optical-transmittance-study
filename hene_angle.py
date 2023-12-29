import numpy as np
import matplotlib.pyplot as plt
from utils.analysis import Measurements, Angles, Maximums, Transmittance
from utils.utils import Index


'''
code to find the best value for the prim's angle from He-Ne data.

we will proceed as follows: first we find the value of alpha that
adjusts the experimental angle associated to the first maximum of T to suit the 
theoretical angle for the same maximum. Next we will do the same for the second one. 
The optimal value of alpha we will choose it to be the mean of both.
'''

# read experimental data

hene_dict = Measurements('data/sources')._hene_tnr()

theta_ext, T_exp = hene_dict['T']

'''
the values to delimiter tha maximums will vary significantly with the prism's angle

the way to obtain the values of the internal angle for the maximums will be to 
search for the maximum values of T in each graph. Since those values won't vary
with the prism's angle, we will find the angle associated to those values by 
finding the index of the transmittance array.

interpolation seems to work better when you need maximums values of T. parabolic 
adjustement doesn't return the right value.
'''

#finding theoretical and experimental values of the transmittance maximums

wavelength = 0.633 
distance = 3.487087087087087e-06 #computed in distance.py
alpha = 45 * np.pi/180 #arbitrary prism's angle to find maximums

theta_int = Angles(alpha = alpha, wavelength = wavelength, 
                   theta_ext = theta_ext).int_angle() # [rad]


# experimental maximums 

plt.figure()
plt.plot((theta_int * 180/np.pi) ,T_exp,'.')
plt.grid(alpha=0.7)
plt.title(r'Experimental transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi))
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()

value_1 = 40.0 ; value_2 = 40.7 ; value_3 = 41.2 #[°] 

T_exp_1 , T_exp_2 , _ , _ = Maximums( T = T_exp ).interpolation(x = (theta_int*180/np.pi), 
                                                           value_1 = value_1, value_2 = value_2, value_3 = value_3)

print('Maximum values of T for experimental graph:', T_exp_1, T_exp_2)

# theoretical maximums

T_teo = Transmittance(distance, alpha, wavelength, theta_ext).transmittance()

plt.figure()
plt.plot((theta_int*180/np.pi), T_teo, '.')
plt.grid(alpha = 0.7)
plt.title(r'Theoretical transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()

value_1 = 40.0 ; value_2 = 40.9 ; value_3 = 41.3 #[°] 

T_teo_1 , T_teo_2 , _ , _ = Maximums( T = T_teo ).interpolation(x = (theta_int*180/np.pi), 
                                                           value_1 = value_1, value_2 = value_2, value_3 = value_3)

print('Maximum values of T for theoretical graph:', T_teo_1, T_teo_2)

'''
the following loop finds the values of alpha that minimizes the distance between 
the theoretical angle where the maximum is and the experimental angle for the 
same maximum
'''

a = np.linspace(44.8, 45.3, 1000) * np.pi/180

min1 = 1000
min2 = 1000

for alpha in a:

    theta_int = Angles(alpha = alpha, wavelength = wavelength, 
                   theta_ext = theta_ext).int_angle() # [rad]
    
    i_1_exp, _ = Index(T_exp).get_double_index(T_exp_1) #the inferior maximum is bivaluated
    i_2_exp = Index(T_exp).get_index(T_exp_2)

    angle_1_exp = (theta_int[i_1_exp]* 180/np.pi) #highest value of theta_int
    angle_2_exp = (theta_int[i_2_exp]* 180/np.pi)

    # plotting the graph to check that everything is right

    # plt.figure()
    # plt.plot((theta_int*180/np.pi),T_exp,'.')
    # plt.grid(alpha=0.7)
    # plt.title(r'Experimental transmittance for $\alpha =  %1.3f°$.' % (alpha*180/np.pi))
    # plt.xlabel(r'$\theta_{int} (°)$')
    # plt.ylabel('T')
    # plt.show()

    T_teo = Transmittance(distance, alpha, wavelength, theta_ext).transmittance()

    i_1_teo, i_2_teo = Index(T_teo).get_double_index(T_teo_1) #theoretical maximums are both equal
    
    angle_1_teo = (theta_int[i_1_teo]* 180/np.pi) #highest value of theta_int
    angle_2_teo = (theta_int[i_2_teo]* 180/np.pi)


    # plt.figure()
    # plt.plot((theta_int*180/np.pi),T_teo,'.')
    # plt.grid(alpha=0.7)
    # plt.title(r'Theoretical transmittance for $\alpha = %1.3f°$' % (alpha*180/np.pi))
    # plt.xlabel(r'$\theta_{int} (°)$')
    # plt.ylabel('T')
    # plt.show()


    if abs(angle_1_exp-angle_1_teo) < min1: #there may be more than one value that minimize the distance between the points ? 
        final_a_1 = alpha
        min1 = abs(angle_1_exp-angle_1_teo)
    else:
        pass

    if abs(angle_2_exp-angle_2_teo) < min2:
        final_a_2 = alpha
        min2 = abs(angle_2_exp-angle_2_teo)
    else:
        pass

print(final_a_2*180/np.pi)
print(final_a_2*180/np.pi)
