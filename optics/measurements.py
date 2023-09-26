import numpy as np
import os


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

