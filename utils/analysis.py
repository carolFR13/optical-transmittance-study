import numpy as np
import os
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from .utils import n_air,n_glass, parabola, Index
from scipy import constants as sc
from uncertainties import umath

class Angles:
    '''
    Angles class to obtain internal angles from the external ones 
    and viceversa.

    :param alpha: prism's angle [rad]
    :param waveleght: light's wavelenght [um]
    :param theta_ext: external angle (optional) [rad]
    :param theta_int: internal angle (optional) [rad]
    
    '''

    def __init__(self, alpha: float, wavelength: float | list, theta_ext: float| list | None = None, 
                 theta_int: float | list | None = None) -> None:
        self.alpha = alpha 
        self.wavelength = wavelength
        self.theta_ext = theta_ext
        self.theta_int = theta_int

        self.n_air = n_air(self.wavelength)
        self.n_glass = n_glass(self.wavelength)

        return None 

    def int_angle(self, uncertainty: bool = False):
        '''
        method to compute the internal angle from the external one
        (the variable theta_ext needs to be defined)

        it allows you to compute internal angle with its uncertainty given an external 
        angle in ufloat format
        '''
        try:
            if uncertainty:
                theta_1 = umath.asin((self.n_air *umath.sin(self.theta_ext))/self.n_glass)
                theta_2 = self.alpha - theta_1
                return theta_2

            else:
                theta_1 = np.arcsin((self.n_air *np.sin(self.theta_ext))/self.n_glass)
                theta_2 = self.alpha - theta_1
                return theta_2
        
        except Exception:
            print('The variable theta_ext is not defined.')
            return None
    
    def int_critical_angle(self):
        '''
        method to compute the internal critical angle
        '''

        theta_c = np.arcsin(self.n_air/self.n_glass)

        return theta_c

    def ext_angle(self):
        '''
        method to compute the external angle from the internal one
        (the variable theta_int needs to be defined)
        '''

        try:
            theta_1=self.alpha - self.theta_int
            theta_ext = np.arcsin((self.n_glass * np.sin(theta_1))/self.n_air)
        except Exception: 
            print('The variable theta_int is not defined.')

        return theta_ext
    
    def critical_lambda(self):
        '''
        if correct method to obtain the critical wavelength
        '''

        def f(wavelength):

            theta_int =  Angles(self.alpha, wavelength, theta_ext = self.theta_ext).int_angle()
            obj = Angles(self.alpha, wavelength, self.theta_ext, theta_int)
            
            theta_2 = obj.int_angle()
            theta_c = obj.int_critical_angle()

            if theta_2 is None: 
                return float('inf')
            else:
                return theta_2 - theta_c

        lambda_c = fsolve(f,x0 = 0.3)

        return lambda_c[0]

class Transmittance: 
    '''
    class to compute the theoretical transmittance and reflection.

    :param distance: distance between prisms [m]
    :param alpha: prism's angle [rad]
    :param wavelength: light's wavelength [um]
    :param theta_ext: external angle [rad]
    '''

    def __init__(self, distance, alpha, wavelength, theta_ext):

        self.distance = distance
        self.alpha = alpha
        self.wavelength = wavelength
        self.theta_ext = theta_ext 

        self.theta_int = Angles(self.alpha, self.wavelength, self.theta_ext).int_angle()
        self.theta_c = Angles(self.alpha, self.wavelength, self.theta_ext).int_critical_angle()

        self.n_air = n_air(self.wavelength)
        self.n_glass = n_glass(self.wavelength)

        return None

    def _fabry_perot(self,_distance, _wavelength, _theta_int, _n_air, _n_glass):
        '''
        internal method to compute trasmission outside TFIR range.
        '''
    
        w = (2*np.pi*sc.c)/(_wavelength*10**(-6))
        beta = _n_glass * w * np.sin(_theta_int)/sc.c

        k1=np.sqrt((w**2*_n_glass**2)/sc.c**2-beta**2)
        k2=np.sqrt((w**2*_n_air**2)/sc.c**2-beta**2)

        T = (1+(1/4)*(((k1**2-k2**2)/(k1*k2))**2)*np.sin(k2*_distance)**2)**(-1)
        R = 1 - T

        return T, R

    def _frustrated_total_int_reflection(self,_distance, _wavelength, _theta_int, _n_air, _n_glass):
        '''
        internal method to compute transmission in TFIR range.
        '''

        w = (2*np.pi*sc.c)/(_wavelength*10**(-6))
        beta = _n_glass * w * np.sin(_theta_int)/sc.c

        k1=np.sqrt((w**2*_n_glass**2)/sc.c**2-beta**2)
        kappa=np.sqrt(beta**2-(w**2*_n_air**2)/sc.c**2)

        T = (1+((kappa**2+k1**2)**2/(4*(k1*kappa)**2))*np.sinh(kappa*_distance)**2)**(-1)
        R = 1 - T

        return T, R
    
    def transmittance(self): 
        '''
        method to compute total transmittance either to a given 
        array of angles or to a given array of wavelegths.

        not defined to compute a single point. 
        '''

        T = []

        if isinstance(self.wavelength,(int,float)):
            for i in range(len(self.theta_int)):
                if self.theta_int[i] < self.theta_c:
                    T.append(self._fabry_perot(self.distance,self.wavelength, 
                                               self.theta_int[i],self.n_air,self.n_glass)[0])
                else:
                    T.append(self._frustrated_total_int_reflection(self.distance,self.wavelength, 
                                                                   self.theta_int[i],self.n_air,self.n_glass)[0])

        elif isinstance(self.wavelength, (list, np.ndarray)):
            for i in range(len(self.wavelength)):
                if self.theta_int[i] < self.theta_c[i]:
                    T.append(self._fabry_perot(self.distance,self.wavelength[i], 
                                               self.theta_int[i],self.n_air[i],self.n_glass[i])[0])
                else:
                    # print('lambda:',self.wavelength[i])
                    T.append(self._frustrated_total_int_reflection(self.distance,self.wavelength[i], 
                                                                   self.theta_int[i],self.n_air[i],self.n_glass[i])[0])
        else:
            print("Invalid type for wavelength")
        
        return T
    
    def reflectance(self): 
        '''
        method to compute total reflectance either to a given 
        array of angles or to a given array of wavelegths.

        not defined to compute a single point. 
        '''

        R = []

        if isinstance(self.wavelength,(int,float)):
            for i in range(len(self.theta_int)):
                if self.theta_int[i] < self.theta_c:
                    R.append(self._fabry_perot(self.distance,self.wavelength, 
                                               self.theta_int[i],self.n_air,self.n_glass)[1])
                else:
                    R.append(self._frustrated_total_int_reflection(self.distance,self.wavelength, 
                                                                   self.theta_int[i],self.n_air,self.n_glass)[1])
        elif isinstance(self.wavelength, (list, np.ndarray)):
        
            for i in range(len(self.wavelength)):
                if self.theta_int[i] < self.theta_c[i]:
                    R.append(self._fabry_perot(self.distance,self.wavelength[i], 
                                               self.theta_int[i],self.n_air[i],self.n_glass[i])[1])
                else:
                    R.append(self._frustrated_total_int_reflection(self.distance,self.wavelength[i], 
                                                                   self.theta_int[i],self.n_air[i],self.n_glass[i])[1])
        else:
            print("Invalid type for wavelength")
        
        return R

class Measurements: 
    '''
    class to read meassured data and store it in dictionaries.

    :param folder_path: name of the folder where the data is stored
    '''

    def __init__(self,folder_path):
        self.folder_path = folder_path

        #dictionaries to store meassured data
        self.lambda_transm_data = {} 
        self.lambda_reflect_data = {}

        self.hene_data = {}

        #dictionaries to store transmittance and reflectance
        self.lambda_transmittance = {}
        self.lambda_reflectance = {}

        self.hene_tnr = {}

    def read_lambda(self, file_path):
        '''
        method to read the meassured data for a wide wavelength range.
        '''

        wavelength = [] ; T = []

        with open(file_path, 'r') as file:
            for l in file: 
                line = l.strip()
                if not line.startswith('"') and line:
                    values = line.split()
                    wavelength.append(float(values[0]))
                    T.append(float(values[1]))

        return np.array(wavelength), np.array(T)

    def read_hene(self,file_path):
        '''
        method to read the messured data as a function of the external angle.
        '''

        P1 = [] ; P2 = [] ; P1_P2 = []
        with open(file_path, 'r') as file:
            for line in file:

                line = line.replace('E','e')
                line = line.replace(',','.')
                values = line.split()

                P1.append(float(values[2]))
                P2.append(float(values[3]))
                P1_P2.append(float(values[4]))

        
        return np.array(P1), np.array(P2), np.array(P1_P2)
    
    def _get_data(self):
        '''
        method to store the directly messured data in the dictionaries
        self.lambda_transm_data, self.lambda_reflect_data and self.hene_data.
        '''

        self.lambda_transm_data = {}
        self.lambda_reflect_data = {}
        self.hene_data = {}

        for file in sorted(os.listdir(self.folder_path)):
            filename = os.path.join(self.folder_path, file)

            if filename.endswith('.txt') and ('SR' in filename):

                _label = os.path.splitext(os.path.basename(filename))[0]
                label = _label[2:]
                wavelength, R = self.read_lambda(file_path = filename)
                self.lambda_reflect_data[label] = [wavelength, R]

            elif filename.endswith('.txt') and ('ST' in filename):

                _label = os.path.splitext(os.path.basename(filename))[0]
                label = _label[2:]
                wavelength, T = self.read_lambda(file_path = filename)
                self.lambda_transm_data[label] = [wavelength, T]

            elif filename.endswith('.txt') and ('video' in filename):
                if 'R' in filename:
                    label = 'R'
                    P1, P2, _ = self.read_hene(filename)
                    self.hene_data[label] = [P1, P2]
                if 'T' in filename:
                    label = 'T'
                    P1, P2, _ = self.read_hene(filename)
                    self.hene_data[label] = [P1, P2]
        return None
        
    def _lambda_transmittance(self):
        '''
        from the meassured data we obtain the transmittance and
        reflectance as T = I_T/(I_T+I_R) and R = I_R/(I_R+I_T)
        being I_R and I_T the messured data in scope mode.
        '''

        self._get_data()

        self.lambda_reflectance = {} ; self.lambda_transmittance = {}

        for label, data in self.lambda_transm_data.items():
            _lambda, T_abs = data
            R_abs = self.lambda_reflect_data[label][1]

            T = T_abs/(R_abs+T_abs)
            R = R_abs/(R_abs+T_abs)

            angle = round((360 - float(label)),1)

            new_label = str(angle)

            self.lambda_transmittance[new_label] = [_lambda, T]
            self.lambda_reflectance[new_label] = [_lambda, R]


        return self.lambda_transmittance
    
    def _lambda_reflectance(self):

        T = self._lambda_transmittance()

        return self.lambda_reflectance
    
    def _hene_tnr(self):
        '''
        method to obtain the transmission and reflection
        from meassured data. In this case P1 is the power of the 
        transmited/reflected light in the prisms and P2 is the reference. 
        Background is the meassured noise.

        the ranges selected to define I_R and I_T were adjusted to 
        correctly fit the external angles defined. 

        returns a dictionary with two tuples R and T with external angle in rad.
        '''

        self._get_data()

        self.hene_tnr = {}

        P1_R, P2_R = self.hene_data['R']
        P1_T, P2_T  = self.hene_data['T']

        background = np.mean(P1_T[1:80])

        I_R = (P1_R[94:613] - background)/P2_R[94:613]
        I_T = (P1_T[100:619] - background)/P2_T[100:619]

        R = I_R /(I_R + I_T)
        T = I_T / (I_R + I_T)

        angles = (360 - np.linspace(356.3,352,len(T))) * np.pi/180

        self.hene_tnr['R'] = [angles, R]
        self.hene_tnr['T'] = [angles, T]

        return self.hene_tnr

class Maximums: 

    '''
    Method to compute the maximas of T and R and the corresponding angles/wavelength 

    :param T: transmittance
    :param R: reflectance
    '''
    def __init__(self, T : None | list = None, R: None | list = None ) -> None:
        self.T = T
        self.R = R

        self.x_R_1: float | None = None #either the value of lambda or theta_int for the 1st max
        self.x_R_2: float | None = None

        self.x_T_1: float | None = None
        self.x_T_2: float | None = None

        self.T_1: float | None = None #value of T for the 1st max
        self.T_2: float | None = None

        self.R_1: float | None = None
        self.R_2: float | None = None
  
    def interpolation(self, x, value_1, value_2, value_3):
        '''
        Method to obtain maximums values of T,R and either the wavelenghts or the
        angles (denoted by x) associated to that maximums by doing a quadratic 
        interpolation.

        value_1, value_2, value_3 are the x values that delimiter the two maximums
        we are studying. x can be either an angles array or a wavelength array. 

        both of them are continuosly increasing functions so we can use get_index()

        recommended to find the maximums of the transmittance.
        '''

        i_1 = Index(x).get_index(value_1)
        i_2 = Index(x).get_index(value_2)
        i_3 = Index(x).get_index(value_3)

        if i_1 > i_3:
            i_1_0 = i_1
            i_1 = i_3
            i_3 = i_1_0

        vect = [self.T]

        if self.R is not None:
            vect.append(self.R)

        for vector in vect:

            vector_1 = np.array(vector[i_1:i_2])
            vector_2 = np.array(vector[i_2:i_3])

            if 'R' in vector:
                index_1 = np.where(vector == vector_1.min())[0][0]
                index_2 = np.where(vector == vector_2.min())[0][0]
            else:
                index_1 = np.where(vector == vector_1.max())[0][0]
                index_2 = np.where(vector == vector_2.max())[0][0]

            v_1 = vector[index_1-10:index_1+10] ; x_1 = x[index_1-10:index_1+10]
            v_2 = vector[index_2-10:index_2+10] ; x_2 = x[index_2-10:index_2+10]
            
            interp_parabolic_1 = interp1d(x_1, v_1, kind='quadratic')
            interp_parabolic_2 = interp1d(x_2, v_2, kind='quadratic')

            #first maximum: 
            x_int_1 = np.linspace(x[index_1-10],x[index_1+9], 10000)
            v_int_1 = interp_parabolic_1(x_int_1)

            #second maximum:
            x_int_2 = np.linspace(x[index_2-10],x[index_2+9], 10000)
            v_int_2 = interp_parabolic_2(x_int_2)

            #store data:

            if 'R' in vector:
                self.R_1 = v_int_1.min()
                self.R_2 = v_int_2.min()
                self.x_R_1 = x_int_1[np.where(v_int_1 ==  v_int_1.min())[0][0]] 
                self.x_R_2 = x_int_2[np.where(v_int_2 ==  v_int_2.min())[0][0]] 
            else:
                self.T_1 = v_int_1.max()
                self.T_2 = v_int_2.max()
                self.x_T_1 = x_int_1[np.where(v_int_1 ==  v_int_1.max())[0][0]] 
                self.x_T_2 = x_int_2[np.where(v_int_2 ==  v_int_2.max())[0][0]] 
            
        return self.T_1, self.T_2, self.x_T_1, self.x_T_2

    def parabolic(self, x, value_1, value_2, value_3):
        '''
        Method to obtain maximums values of T,R and either the wavelenghts or the
        angles (denoted by x) associated to that maximums by doing a parabolic adjust.

        value_1, value_2, value_3 are the x values that delimiter the two maximums
        we are studying. x can be either an angles array or a wavelength array. 

        both of them are continuosly increasing functions so we can use get_index()
        '''

        i_1 = Index(x).get_index(value_1)
        i_2 = Index(x).get_index(value_2)
        i_3 = Index(x).get_index(value_3)

        if i_1 > i_3:
            i_1_0 = i_1
            i_1 = i_3
            i_3 = i_1_0

        vect = [self.T]
        if self.R is not None:
            vect.append(self.R)

        for vector in vect:

            vector_1 = np.array(vector[i_1:i_2])
            vector_2 = np.array(vector[i_2:i_3])

            if 'R' in vector:
                index_1 = np.where(vector == vector_1.min())[0][0]
                index_2 = np.where(vector == vector_2.min())[0][0]
            else:
                index_1 = np.where(vector == vector_1.max())[0][0]
                index_2 = np.where(vector == vector_2.max())[0][0]

            v_1 = vector[index_1-10:index_1+10] ; x_1 = x[index_1-10:index_1+10]
            v_2 = vector[index_2-10:index_2+10] ; x_2 = x[index_2-10:index_2+10]
            
            params1, _ = curve_fit(parabola, x_1, v_1)
            params2, _ = curve_fit(parabola, x_2, v_2)

            #first maximum: 
            x_int_1 = np.linspace(x[index_1-10],x[index_1+9], 10000)
            v_int_1 = parabola(x_int_1, *params1)

            #second maximum:
            x_int_2 = np.linspace(x[index_2-10],x[index_2+9], 10000)
            v_int_2 = parabola(x_int_2, *params2)

            #store data:

            if 'R' in vector:
                self.R_1 = v_int_1.min()
                self.R_2 = v_int_2.min()
                self.x_R_1 = x_int_1[np.where(v_int_1 ==  v_int_1.min())[0][0]] 
                self.x_R_2 = x_int_2[np.where(v_int_2 ==  v_int_2.min())[0][0]] 
            else:
                self.T_1 = v_int_1.max()
                self.T_2 = v_int_2.max()
                self.x_T_1 = x_int_1[np.where(v_int_1 ==  v_int_1.max())[0][0]] 
                self.x_T_2 = x_int_2[np.where(v_int_2 ==  v_int_2.max())[0][0]] 
            
        return self.T_1, self.T_2, self.x_T_1, self.x_T_2
        
class Angles_max(Maximums):
    '''
    class to obtain prism's distance when trasnmitance is a function of 
    theta_int
    '''

    def __init__(self, wavelength: float, T : None | list = None, R: None | list = None ) -> None:

        super().__init__(T = T, R = R)

        self.wavelength = wavelength
    
    def distance(self, theta_1: float | None = None, theta_2: float | None = None,
                 theta: list | None = None, value_1: float | None = None, 
                 value_2: float | None = None, value_3: float | None = None ) -> float:
        
        '''
        method to compute the distance. You either give the 2 values of the angles 
        corresponding to the maximums in T/R (theta_1, theta_2) or either
        give the angles array (theta) and the values value_1, value_2, value_3
        to compute the 2 angles corresponding to the maxima.

        '''

        if (theta_1, theta_2) == (None, None):

            try:
                _ , _ , theta_1, theta_2 = self.interpolation(x = theta, value_1 = value_1, 
                                                            value_2 = value_2, value_3 = value_3)
            except Exception:

                print("You didn't provide the right arguments" )
            
        n_air_ = n_air(self.wavelength)
        n_glass_ = n_glass(self.wavelength)

        d = (1/2) * ( (( n_air_**2 - (n_glass_*np.sin(theta_2))**2 )**(1/2))/self.wavelength - (( n_air_**2 - (n_glass_*np.sin(theta_1))**2 )**(1/2))/self.wavelength  )**(-1) 
            
        return abs(d) 

class Lambda_max(Maximums):
    '''
    class to obtain prism's distance when transmittance is a wavelength's 
    function
    '''

    def __init__(self, wavelength: list, T : None | list = None, R: None | list = None ) -> None:

        super().__init__(T,R)

        self.wavelength = wavelength

    def distance(self, theta_1: float | None = None, theta_2: float | None = None,
                 lambda_1: float | None = None, lambda_2: float | None = None,
                 alpha: float | None = None, theta_ext: float | None = None,
                 value_1: float | None = None, value_2: float | None = None, 
                 value_3: float | None = None, uncertainty: bool = False) -> float:
        
        '''
        method to compute the distance. you either provide the internal angles and wavelengths 
        values for the T maximums (theta_1, theta_2, lambda_1, lambda_2) 
        or either the values of the prism's angle, theta_ext and values 1,2,3 to compute the internal 
        angles and wavelenghts corresponding to the maxima.
        '''

        if (theta_1, theta_2, lambda_1, lambda_2) == (None, None, None, None):

            try:
                _ , _ , lambda_1, lambda_2 = self.interpolation(x = self.wavelength, value_1 = value_1, 
                                                            value_2 = value_2, value_3 = value_3)
                
                theta_1 = Angles(alpha = alpha, wavelength= lambda_1, theta_ext = theta_ext).int_angle(uncertainty = uncertainty)
                theta_2 = Angles(alpha = alpha, wavelength= lambda_2, theta_ext = theta_ext).int_angle(uncertainty = uncertainty)

            except Exception:
                print("You didn't provide the right arguments" )

        n_air_1 = n_air(lambda_1) ; n_air_2 = n_air(lambda_2)
        n_glass_1 = n_glass(lambda_1) ; n_glass_2 = n_glass(lambda_2)

        if uncertainty:
            d = (1/2) * ( (( n_air_2**2 - (n_glass_2*umath.sin(theta_2))**2 )**(1/2))/lambda_2 - (( n_air_1**2 - (n_glass_1*umath.sin(theta_1))**2 )**(1/2))/lambda_1  )**(-1) 
        
        else:
            d = (1/2) * ( (( n_air_2**2 - (n_glass_2*np.sin(theta_2))**2 )**(1/2))/lambda_2 - (( n_air_1**2 - (n_glass_1*np.sin(theta_1))**2 )**(1/2))/lambda_1  )**(-1) 

        return abs(d)
