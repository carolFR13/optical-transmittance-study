import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



'''
In this file we define the classes to read the data from 
.txt files in order to store it in a dictionary.
'''

class Wide_range:
    '''
    class to read measured data for a wide spectrum.

    :param folder_path: name of the folder where the data is stored
    '''

    def __init__(self, folder_path):

        self.folder_path = folder_path

        # dictionary to store measured data
        self.measured_transmision = {}
        self.measured_reflection = {}

        #dictionaries to store transmittance and reflectance
        self.transmittance = {}
        self.reflectance = {}
        
    

    def read(self, file_path):
        '''
        method to read the messured data for a wide spectrum
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
    

    def _get_data(self):
        '''
        method to store the directly messured data in the dictionaries.

        '''

        self.measured_transmision = {}
        self.measured_reflection = {}

        for file in sorted(os.listdir(self.folder_path)):
            filename = os.path.join(self.folder_path, file)

            if filename.endswith('.txt') and ('SR' in filename):

                _label = os.path.splitext(os.path.basename(filename))[0]
                label = _label[2:]
                wavelength, R = self.read(file_path = filename)
                self.measured_reflection[label] = [wavelength, R]

            elif filename.endswith('.txt') and ('ST' in filename):

                _label = os.path.splitext(os.path.basename(filename))[0]
                label = _label[2:]
                wavelength, T = self.read(file_path = filename)
                self.measured_transmision[label] = [wavelength, T]
        
        return None
    
    def transmittance_computation(self):
        '''
        from the measured data we obtain the transmittance and
        reflectance as T = P_T/(P_T+P_R) and R = P_R/(P_R+P_T)
        being P_R and P_T the messured data in scope mode.
        '''

        self._get_data()

        self.reflectance = {} ; self.transmittance = {}

        for label, data in self.measured_transmision.items():
            _lambda, T_abs = data
            R_abs = self.measured_reflection[label][1]

            T = T_abs/(R_abs+T_abs)
            R = R_abs/(R_abs+T_abs)

            angle = round((360 - float(label)),1)

            new_label = str(angle)

            self.transmittance[new_label] = [_lambda, T]
            self.reflectance[new_label] = [_lambda, R]

        return self.transmittance
    
    def reflectance_computantion(self):

        T = self.transmittance_computation()

        return self.reflectance


class Monochromatic:

    '''
    Class to read measured data from a monochromatic laser. 
    
    There are two ways of measuring the transmitted and reflected power:

    1. Using a glass to measure a reference, so we obtain either the relectance
    or the transmittance by computing the fraction: trasnmitted (reflected) power /
    reference.

    2. measuring simultaneously the transmitted and reflected power so we obtain
    in the same file R and T.

    '''

    def __init__(self, folder_path):
        
        self.folder_path = folder_path

        #  dictionary to store measured P1 and P2
        self.measured_data = {}

        # dictionary to store values to adjust
        self.values = {}

        # dictionary to store detector constants for 2nd method
        self.constants = {}

        # dictionary to store corrected values for 2nd method
        self.corrected_data = {}

        # dictionary to store T and R computed from directly measurements.
        self.data = {}
        

    def read(self, file_path):
        '''
        It returns the arrays with the measured values without computing anything.

        '''

        P1 = [] ; P2 = []
        with open(file_path, 'r') as file:
            for line in file:

                line = line.replace('E','e')
                line = line.replace(',','.')
                values = line.split()

                P1.append(float(values[2]))
                P2.append(float(values[3]))

        return np.array(P1), np.array(P2)
    
    def _get_data(self, reference : bool = False):

        '''
        Method to store the data from the files as it was measured in a
        dictionary.

        The data must be stored in a folder where we have 4 files 
        in the form:

        1. For the first method we will have the files TE_R, TE_T, 
        TM_R, TM_T

        2. For the second method the files would be TE_1, TE_2, TM_1, 
        TM_2. In this case we use the format explained in 'name_format.txt'.

        '''

        # if you want your dictionary not to accumulate the data 
        # from the different files uncomment the next line.

        self.measured_data = {}

        if reference:

            for file in sorted(os.listdir(self.folder_path)):

                filename = os.path.join(self.folder_path, file)

                if filename.endswith('.txt') and ('TE' or 'TM' in filename):

                    basename = os.path.basename(filename)

                    label_1 = (os.path.splitext(basename)[0]).split('_')[0]
                    label_2 = (os.path.splitext(basename)[0]).split('_')[1]

                    P1, P2 = self.read(filename)

                    self.measured_data[label_1][label_2] = [P1, P2]
        
        else: 
                '''
                We will store the measured data in a dictionary. According
                to our name format the files are organised with a 1 or 2 label
                that corresponds to the order of the detectors as follows:

                    1 : T -> S120UV ; R -> S121B
                    2 : T -> S121B ; R -> S120UV

                Here we will change those labels (1,2) in order to adjust them
                to correspond with the detector used as follows:

                    1 : detector S120UV (associated to channel 1, that is P1)
                    2 : detector S121B (associated to channel 2, that is P2)
                
                '''


                for file in sorted(os.listdir(self.folder_path)):

                    filename = os.path.join(self.folder_path, file)

                    if filename.endswith('.txt') and ('data_info' not in filename):

                        basename = os.path.basename(filename)
                        labels = (os.path.splitext(basename)[0]).split('_')

                        # 'a' if medium is air ; 'b' if medium is water
                        label_1 = labels[0] 

                        # 'a' if thickness provided by Alejandro ; 'r' if thickness provided by Raúl
                        label_2 = labels[1]

                        # TE ; TM
                        label_3 =  labels[2]

                        # 1 ; 2
                        label_4 = labels[3]

                        _label = label_1 + '_' + label_2 + '_' + label_3

                        if label_4 == '1':
                            T1, R2 = self.read(filename)
                            if _label not in self.measured_data:
                                self.measured_data[_label] = {'T': {}, 'R': {}}
                            self.measured_data[_label]['T']['1'] = T1
                            self.measured_data[_label]['R']['2'] = R2
                        else:
                            R1, T2 = self.read(filename)
                            if _label not in self.measured_data:
                                self.measured_data[_label] = {'T': {}, 'R': {}}
                            self.measured_data[_label]['T']['2'] = T2
                            self.measured_data[_label]['R']['1'] = R1

        return None
    
    def _adjust_data(self):

        folderpath = self.folder_path
        filepath = os.path.join(folderpath, 'data_info.txt')

        self._get_data()

        with open(filepath, 'r') as file:
            data = file.readlines()

        self.values = {}
        for i in range(2,len(data),4):
            line = data[i].strip()
            if line:
                file = line.split('\t')[1]
                labels  = file.split('_')
                label =  labels[0] + '_' + labels[1] + '_' + labels[2]
                label2 = labels[3]

                values = data[i+1].replace('E','e').replace(',','.')
                values_start = values.strip().split('\t')

                values2 = data[i+2].replace('E','e').replace(',','.')
                values_end = values2.strip().split('\t')


                if label not in self.values:
                    self.values[label] = {'T': {}, 'R': {}}
                
                if label2 == '1':
                    self.values[label]['T']['1'] = (values_start[1], values_start[3], values_end[1], values_end[3])
                    self.values[label]['R']['2'] = (values_start[2], values_start[3], values_end[2], values_end[3])
                elif label2 == '2':
                    self.values[label]['R']['1'] = (values_start[1], values_start[3], values_end[1], values_end[3])
                    self.values[label]['T']['2'] = (values_start[2], values_start[3], values_end[2], values_end[3])
                else:
                    print('Something may be wrong')

        for key, dict in self.measured_data.items():
            for key_2, dict_2 in dict.items():
                for key_3 in ('1','2'):

                    vector = list(self.measured_data[key][key_2][key_3])
                    init_value, init_angle, end_value, end_angle = self.values[key][key_2][key_3]

                    index_1 = vector.index(float(init_value))
                    index_2 = vector.index(float(end_value))
                    vector_n = vector[index_1:index_2]
                    print(key, key_2, key_3, '      ',len(vector_n))

                    angles = np.linspace(float(init_angle),float(end_angle),len(vector_n))

                    self.measured_data[key][key_2][key_3] = (vector_n, angles)

        return None

    
    # Specific methods for the 2nd way of measurement

    def _get_constants(self):

        '''
        Method to obtain the constants associated to each detector 
        if different detectors are used to measure simultaneously R and T.

        From measured power we will get the fraction k1/k2, 
        being k1,2 the constant associated to the 1,2 detector.
        '''

        self._adjust_data()

        # self.constants = {}

        try:

            for key_1, _ in self.measured_data.items():

                T_1 = self.measured_data[key_1]['T']['1'][0]
                T_2 = self.measured_data[key_1]['T']['2'][0]
                R_1 = self.measured_data[key_1]['R']['1'][0]
                R_2 = self.measured_data[key_1]['R']['2'][0]

                k1_k2 = np.sqrt((np.multiply(R_1,T_1))/np.multiply(R_2,T_2))

                self.constants[key_1] = k1_k2 

        except ValueError:
            print('Meassurements do not have the same number of elements for', str(key_1))

        return None
    
    def _corrected_values(self):

        self._adjust_data()
        self._get_constants()

        self.corrected_data = {}

        for key_1, _  in self.measured_data.items():

            T_1 , a_1 = self.measured_data[key_1]['T']['1']
            T_2 , a_2 = self.measured_data[key_1]['T']['2']
            R_1 , b_1 = self.measured_data[key_1]['R']['1']
            R_2 , b_2 = self.measured_data[key_1]['R']['2']

            k1_k2 = self.constants[key_1]

            T_1_n = T_1 * (1/k1_k2)
            T_2_n = T_2 * k1_k2
            R_1_n = R_1 * (1/k1_k2)
            R_2_n = R_2 * k1_k2

            if key_1 not in self.corrected_data:
                self.corrected_data[key_1] = {'T': {}, 'R': {}}

            self.corrected_data[key_1]['T']['1'] = T_1_n , a_1
            self.corrected_data[key_1]['R']['2'] = R_2_n , b_2
            self.corrected_data[key_1]['T']['2'] = T_2_n , a_2
            self.corrected_data[key_1]['R']['1'] = R_1_n , b_1

        return None


    def _get_transmittance(self, reference : bool = False):

        self._get_data(reference = reference)

        # self.data={}

        if reference: 

            for  key_1, _ in self.measured_data.items():

                P1_R , P2_R = self.measured_data[key_1]['R']
                P1_T , P2_T = self.measured_data[key_1]['T']

                P_R = P1_R / P2_R
                P_T = P1_T / P2_T

                R = P_R / (P_R + P_T)
                T = P_T / (P_R + P_T)
            
                self.data[key_1]['R'] = R
                self.data[key_1]['T'] = T

        else: 

            '''
            Here we are going to compute the trasnmittance and
            reflectance from the measured data. We have aready
            defined a method to correct the data in order to not 
            be dependent of the detector used. 

            As we have defined previously our data is organiced as:

                1 : detector S120UV (associated to channel 1, that is P1)
                2 : detector S121B (associated to channel 2, that is P2)
                
            Since now we want to obtain T and R as T = P_T / (P_R + P_T)
            and R = P_R / (P_T + P_R) in order to eliminate the fluctuations
            of the laser, we need to go back to the previous notation where:

                1 : T -> S120UV ; R -> S121B
                2 : T -> S121B ; R -> S120UV

            Defined like this we will have a T and R associated to the 
            measurement and not to the detector. The computantion we will do 
            would be: 

                1 : T = T_1 / (T_1 + R_2) and R = R_2 / (T_1 + R_2)
                2 : T = T_2 / (T_2 + R_1) and R = R_1 / (T_2 + R_1)

            and if everything is correct we will have the same values for 
            both measurements.

            '''

            self._get_constants()
            self._corrected_values()

            for key_1, _  in self.corrected_data.items():

                T_1 , a_1 = self.corrected_data[key_1]['T']['1']
                T_2 , a_2 = self.corrected_data[key_1]['T']['2']
                R_1 , b_1 = self.corrected_data[key_1]['R']['1']
                R_2 , b_2 = self.corrected_data[key_1]['R']['2']

                T_1_n = T_1 / (T_1 + R_2)
                R_1_n = R_2 / (T_1 + R_2)

                T_2_n = T_2 / (T_2 + R_1)
                R_2_n = R_1 / (T_2 + R_1)

                if key_1 not in self.data:
                    self.data[key_1] = {'T': {}, 'R': {}}

                self.data[key_1]['T']['1'] = T_1_n , a_1
                self.data[key_1]['R']['2'] = R_2_n , b_1
                self.data[key_1]['T']['2'] = T_2_n , a_2
                self.data[key_1]['R']['1'] = R_1_n , b_2

        return None




def report(folder_path, file, limit):

    obj = Monochromatic(folder_path)
    obj._adjust_data()
    obj._get_constants()
    obj._corrected_values()
    obj._get_transmittance()

    with PdfPages(file) as pdf:

        for key, _ in obj.measured_data.items():
            for key2 in ['T', 'R']:
                matrix = []
                K1_K2 = np.array(obj.constants[key][:limit])
                P_1 = np.array(obj.measured_data[key][key2]['1'][0][:limit])
                P_1_n = np.array(obj.corrected_data[key][key2]['1'][0][:limit])
                P_2 = np.array(obj.measured_data[key][key2]['2'][0][:limit])
                P_2_n = np.array(obj.corrected_data[key][key2]['2'][0][:limit])
                if key2 == 'T':
                    T_1 = np.array(obj.data[key][key2]['1'][0][:limit])
                    T_2 = np.array(obj.data[key][key2]['2'][0][:limit])
                    col_head = [r'$K_1/K_2$',r'$P_{T,1}$',r'$P^c_{T,1}$', r'$T_1$',r'$P_{T,2}$', r'$P^c_{T,2}$', r'$T_2$']

                else:
                    T_1 = np.array(obj.data[key][key2]['2'][0][:limit])
                    T_2 = np.array(obj.data[key][key2]['1'][0][:limit])
                    col_head = [r'$K_1/K_2$',r'$P_{R,1}$',r'$P^c_{R,1}$', r'$R_1$',r'$P_{R,2}$',r'$P^c_{R,2}$', r'$R_2$']

                
                matrix.append([K1_K2, P_1, P_1_n, T_1, P_2, P_2_n, T_2])
                matrix = np.array(matrix[0])
            
                table = [[f"{num:.3e}" for num in sublist] for sublist in matrix.T]


                fig, ax = plt.subplots(figsize=(8,6))
                ccolors = np.full(len(col_head), 'lavender')

                tb = plt.table(cellText = table,
                                    colLabels=col_head,
                                    rowLabels=None,
                                    colColours=ccolors,
                                    cellLoc='center',
                                    loc='center')
                tb.scale(1, 2)
                tb.set_fontsize(16)
                ax.axis('off')
                plt.title(key+'_'+key2)
                pdf.savefig(fig)
                plt.close()
            

            # Reflected power plot
            plt.figure()
            P_R_1 , angle_1 = obj.measured_data[key]['R']['1']
            P_R_2 , angle_2 = obj.measured_data[key]['R']['2']
            P_R_1_c , angle_1 = obj.corrected_data[key]['R']['1']
            P_R_2_c , angle_2 = obj.corrected_data[key]['R']['2']
            plt.title('Reflected power for '+key+' file')
            plt.grid(alpha=0.7)
            plt.plot(angle_1,P_R_1,label=r'$P_{R,1}$')
            plt.plot(angle_2,P_R_2,label=r'$P_{R,2}$')
            plt.plot(angle_1,P_R_1_c,label=r'$P^c_{R,1}$')
            plt.plot(angle_2,P_R_2_c,label=r'$P^c_{R,2}$')
            plt.legend(loc='best')
            plt.xlabel(r'$\theta$(°)')
            plt.ylabel(r'$P_R$')
            pdf.savefig()
            plt.close()

            # Transmitted power plot
            plt.figure()
            P_T_1 , angle_1 = obj.measured_data[key]['T']['1']
            P_T_2 , angle_2 = obj.measured_data[key]['T']['2']
            P_T_1_c , angle_1 = obj.corrected_data[key]['T']['1']
            P_T_2_c , angle_2 = obj.corrected_data[key]['T']['2']
            plt.title('Trasmitted power for '+key+' file')
            plt.grid(alpha=0.7)
            plt.plot(angle_1,P_T_1,label=r'$P_{T,1}$')
            plt.plot(angle_2,P_T_2,label=r'$P_{T,2}$')
            plt.plot(angle_1,P_T_1_c,label=r'$P^c_{T,1}$')
            plt.plot(angle_2,P_T_2_c,label=r'$P^c_{T,2}$')
            plt.legend(loc='best')
            plt.xlabel(r'$\theta$(°)')
            plt.ylabel(r'$P_T$')
            pdf.savefig()
            plt.close()

            # Reflectance plot
            plt.figure()
            R_1 , angle_1 = obj.data[key]['R']['1']
            R_2 , angle_2 = obj.data[key]['R']['2']
            plt.title('Reflectance for '+key+' file')
            plt.grid(alpha=0.7)
            plt.plot(angle_1,R_1,'.',label='1st measurement')
            plt.plot(angle_2,R_2,'.',label='2nd measurement')
            plt.legend(loc='best')
            plt.xlabel(r'$\theta$(°)')
            plt.ylabel('R')
            pdf.savefig()
            plt.close()

            # Transmitance plot
            plt.figure()
            T_1 , angle_1 = obj.data[key]['T']['1']
            T_2 , angle_2 = obj.data[key]['T']['2']
            plt.title('Transmitance for '+key+' file')
            plt.grid(alpha=0.7)
            plt.plot(angle_1,T_1,'.',label='1st measurement')
            plt.plot(angle_2,T_2,'.',label='2nd measurement')
            plt.legend(loc='best')
            plt.xlabel(r'$\theta$(°)')
            plt.ylabel('T')
            pdf.savefig()
            plt.close()

    return None







class Measurements: 
    '''
    Deprecated class. Used for initial analysis. Wide_range and Monochromatic
    classes are now available.

    :param folder_path: name of the folder where the data is stored
    '''

    def __init__(self,folder_path):
        self.folder_path = folder_path

        #dictionaries to store measured data
        self.lambda_transm_data = {} 
        self.lambda_reflect_data = {}

        self.hene_data = {}

        #dictionaries to store transmittance and reflectance
        self.lambda_transmittance = {}
        self.lambda_reflectance = {}

        self.hene_tnr = {}

    def read_lambda(self, file_path):
        '''
        method to read the measured data for a wide wavelength range.
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

        P1 = [] ; P2 = []
        with open(file_path, 'r') as file:
            for line in file:

                line = line.replace('E','e')
                line = line.replace(',','.')
                values = line.split()

                P1.append(float(values[2]))
                P2.append(float(values[3]))

        return np.array(P1), np.array(P2)
    
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
                    P1, P2 = self.read_hene(filename)
                    self.hene_data[label] = [P1, P2]
                if 'T' in filename:
                    label = 'T'
                    P1, P2 = self.read_hene(filename)
                    self.hene_data[label] = [P1, P2]
        return None
        
    def _lambda_transmittance(self):
        '''
        from the measured data we obtain the transmittance and
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
    
    def _hene_tnr(self, start = False, end = False):
        '''
        method to obtain the transmission and reflection
        from measured data. In this case P1 is the power of the 
        transmited/reflected light in the prisms and P2 is the reference. 
        Background is the measured noise.

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

