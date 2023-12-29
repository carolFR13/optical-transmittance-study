import numpy as np
import matplotlib.pyplot as plt
from utils.analysis import Measurements, Angles, Maximums, Transmittance


'''
code to find the distance between prisms with He-Ne data. 

we use experimental graph to find the relative angular distance between 
the two maximums. the distance would be the one that provides the same
relative angular distance in the theorical graph.

we start assuming a particular value of the prism's angle.
'''

# read experimental data

hene_dict = Measurements('data/sources')._hene_tnr()

theta_ext, T_exp = hene_dict['T']

alpha = 45.214414414414414 * np.pi/180
wavelength = 0.633 #HeNe wavelength

theta_int = Angles(alpha = alpha, wavelength = wavelength, theta_ext = theta_ext).int_angle() # [rad]

plt.figure()
plt.plot((theta_int*180/np.pi), T_exp, '.')
plt.grid(alpha = 0.7)
plt.title(r'Experimental transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()

# the best values to delimiter the maximums for this alpha are: 

value_1 = 40.0 # [°]
value_2 = 40.7 
value_3 = 41.3

# compute the relative angular distance for experimental data

_ , _ , theta_1 , theta_2 = Maximums( T = T_exp).parabolic(x = (theta_int*180/np.pi), 
                                                           value_1 = value_1, value_2 = value_2, value_3 = value_3)


relat_dist_exp = abs(theta_1-theta_2) #[°]


'''
we compute the distance by choosing the value that minimizes the 
difference between the relative angular distance for the experimental data 
and the one we obtain from theoretical data.

the relative distance is computed with the values 1,2,3 found for the experimental
graph since we want both graphs to look alike (?)
'''

d = np.linspace(2.6,4,100) * 10**(-6)  #we expect the distance to be near 3 um

min = 1000
possible_d = []
min_vector= []

for dist in d:

    T_teo = Transmittance(dist, alpha, wavelength, theta_ext).transmittance()
    _ , _ , theta_1 , theta_2 = Maximums( T = T_teo).parabolic(x = (theta_int*180/np.pi), 
                                                           value_1 = value_1, value_2 = value_2, value_3 = value_3)
    
    relat_dist_teo = abs(theta_1-theta_2)
    min_vector.append(abs(relat_dist_exp - relat_dist_teo))
    if abs(relat_dist_exp - relat_dist_teo) < min:
        min = abs(relat_dist_exp - relat_dist_teo)
        final_d = dist
    else:
        continue

print(final_d)
print(min_vector)

# final_d = 3.484284284284284e-06 for 45.13°
# final_d = 3.4856856856856854e-06 for 45.21341341341341°
# final_d = 3.487087087087087e-06 for 45.21391391391391°
# final_d = 3.487087087087087e-06 for 45.214414414414414°

T_teo = Transmittance(final_d, alpha, wavelength, theta_ext).transmittance()

plt.figure()
plt.plot((theta_int*180/np.pi), T_teo, '.')
plt.grid(alpha = 0.7)
plt.title(r'Theoretical transmittance for $\alpha$ = %2.2f°' % (alpha*180/np.pi) )
plt.xlabel(r'$\theta_{int} (°)$')
plt.ylabel('T')
plt.show()



plt.figure()
plt.plot(d, min_vector,'.')
plt.grid(alpha = 0.7)
plt.title(r'Obtention of distance between prisms' )
plt.xlabel(r' d ($\mu m$)')
plt.ylabel(r'|$\Delta d_{exp}$-$\Delta d_{teo}$|')
plt.show()


