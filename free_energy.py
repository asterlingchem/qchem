import math
import sys, os
import argparse

def get_args():
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('filenames',
                           nargs='+',
                           help="Input file names, including '.out' suffix")
    my_parser.add_argument("-c", "--cutoff",
                           action='store',
                           type=float,
                           default=100,
                           help='Frequency cut-off for qRRHO analysis in cm^-1. Default: 100')
    my_parser.add_argument("-t", "--temp",
                           action='store',
                           type=float,
                           default=298.15,
                           help='Temperature for free energies in K. Default: 298.15')
    my_parser.add_argument("-i", "--imag",
                           action='store',
                           type=float,
                           default=0,
                           help='Reverse sign of imaginary frequency if below cut-off, in cm^-1. Default: 0')
    my_parser.add_argument("-v", "--verbose",
                           action='store_true',
                           help='Verbose printing')
    my_parser.add_argument("-q", "--partition",
                           action='store_true',
                           help="Print all partition functions")
    return my_parser.parse_args()

# verbose = True

class QC_freq_to_FreeEnergy:
    args = get_args()
    # T = float(sys.argv[1])              # temperature in SI
    T = float(args.temp)              # temperature in SI
    k = 1.38064852e-23      # k in SI
    kT = 1.38064852e-23*T   # kT in SI
    h = 6.62607004e-34      # planck const in SI
    hbar = h/(2*math.pi)
    c = 299792458           # speed of light in SI
    P = 101325              # pressure in SI
    V = 8.314462618*T/P     # volume of 1 mol ideal gas
    NA = 6.02214076e23
    amu_to_si = 1.660539066e-27     # amu to SI
    boh_to_si = 5.291772109e-11     # bohr radius to SI
    moi_to_si = amu_to_si * (boh_to_si)**2  # moment of inertia unit to SI
    Jtokcal = 2.390057e-4
    kcaltoHa = 0.001593601097421
    i_cutoff = float(args.imag)
    def __init__(self, filename):
        self.lines = None
        self.scf_energy = None
        self.data = {'scf_energy':None, 'vib_freq':[], 'mom_inert':[], 'mol_mass':None, 'sym_num':None, 'imag_freq':[]}
        self.vib_ps = None      # partition sum of one molecule
        self.rot_ps = None
        self.trans_ps = None
        self.total_ps = None
        self.vib_E = None       # energy of 1mol molecule
        self.vib_E_quasi_RRHO = None
        self.trans_E = None
        self.rot_E = None
        self.total_E = None
        self.total_E_quasi_RRHO = None
        self.scf_plus_h_corr = None
        self.scf_plus_q_h_corr = None
        self.vib_S = None       # entropy of 1mol molecule
        self.vib_S_quasi_RRHO = None
        self.trans_S = None
        self.rot_S = None
        self.total_S = None
        self.total_S_quasi_RRHO = None
        self.vib_FreeEnergy = None     # free energy of 1mol molecule
        self.vib_FreeEnergy_quasi_RRHO = None
        self.rot_FreeEnergy = None
        self.trans_FreeEnergy = None
        self.total_FreeEnergy = None
        self.total_FreeEnergy_quasi_RRHO = None
        self.scf_plus_g_corr = None
        self.scf_plus_q_g_corr = None

        f = open(filename, 'r')
        self.lines = f.readlines()
        f.close()

    def print_raw_data(self):
        for line in self.lines:
            print(line)

    def get_scf_energy(self):
       '''
       Return electronic energy from SCF calculation
       '''
       for line in self.lines:
           if "Total energy in the final basis set " in line:
               self.data['scf_energy'] = float(line.split()[-1])
               self.scf_energy = float(line.split()[-1])

    def get_useful_data(self):
        freq_data = []
        for line in self.lines:
            if 'Frequency:' in line:        # read vibrational frequency in cm^-1
                freq_data.extend(line.split()[1:])
            if 'Eigenvalues --' in line:    # read principle moment of inertia in amu*Bohr^2
                self.data['mom_inert'].extend(line.split()[2:])
            if 'Molecular Mass' in line:    # read molecular mass in amu
                self.data['mol_mass'] = float(line.split()[-2])
            if 'Rotational Symmetry Number' in line:    # read rotational symmetry number
                self.data['sym_num'] = int(line.split()[-1])

        imag_above_thresh = 0
        imag_below_thresh = 0
        for i in freq_data:
            if float(i) < 0 and abs(float(i)) >= args.imag:
                print(f"Imaginary mode above threshold |{args.imag}| cm^-1 detected")
                self.data['imag_freq'].append(float(i))
                imag_above_thresh += 1
            elif float(i) < 0 and abs(float(i)) < args.imag:
                print(f"Imaginary mode below threshold |{args.imag}| cm^-1 detected")
                print("Reversing sign of mode below threshold")
                self.data['imag_freq'].append(float(i))
                i = abs(float(i))
                self.data['vib_freq'].append(float(i))
                imag_below_thresh += 1
            elif float(i) > 0:
                self.data['vib_freq'].append(float(i))
        print(f"Number of imaginary modes above threshold |{args.imag}| cm^-1: {imag_above_thresh}")
        print(f"Number of imaginary modes below threshold |{args.imag}| cm^-1: {imag_below_thresh}")
        total_imag_modes = imag_above_thresh + imag_below_thresh
        if total_imag_modes > 1:
            if total_imag_modes == imag_above_thresh:
                print("   ******* NOT a true minimum or TS *******")
                print("               Tread carefully             ")
                print(f"    Imaginary frequencies: {self.data['imag_freq']} cm^-1")
            elif total_imag_modes == imag_below_thresh:
                print("   ******* NOT a true minimum or TS *******")
                print("      All imaginary modes below threshold  ")
                print("              Treating as minimum          ")
                print(f"    Imaginary frequencies: {self.data['imag_freq']} cm^-1")
            elif imag_above_thresh == 1:
                print("   ******* NOT a true minimum or TS *******")
                print("     One imaginary mode above threshold    ")
                print("                Treating as TS             ")
                print(f"    Imaginary frequencies: {self.data['imag_freq']} cm^-1")
        elif total_imag_modes == 1:
            if imag_above_thresh == 1:
                print("         ******* FOUND A TS *******        ")
                print("             One imaginary mode            ")
                print("                Treating as TS             ")
                print(f"    Abs TS frequency: {self.data['imag_freq'][0]} cm^-1")
            elif imag_below_thresh == 1:
                print("   ******* NOT a true minimum or TS *******")
                print("             One imaginary mode            ")
                print("            Treating as minimum            ")
                print(f"    Abs lowest frequency: {self.data['imag_freq'][0]} cm^-1")
        else:
            print("    ********* Found a MINIMUM *********    ")

        for i in range(len(self.data['mom_inert'])):
            self.data['mom_inert'][i] = float(self.data['mom_inert'][i])

    def print_useful_data(self):
        if verbose == True:
            for key, val in self.data.items():
                print(key, ' : ', val)
        # else:
        #     lowest_mode = self.data['vib_freq'][0]
        #     second_lowest_mode = self.data['vib_freq'][1]
        #     print("       Lowest frequency (cm-1):", lowest_mode)
        #     print("Second lowest frequency (cm-1):", second_lowest_mode)
        #     if lowest_mode >= 0:
        #         print("   ********* Found a MINIMUM *********")
        #     elif second_lowest_mode <= 0:
        #         print("   ******* NOT a minimum or TS *******")
        #     else:
        #         print("*** One imaginary mode – TREATING AS TS ***")
    # calculate vibrat1ional partition sum of one molecule
    def cal_vib_ps(self):
        self.vib_ps = 1
        for freq in self.data['vib_freq']:
            ps = math.exp(-self.h*self.c*100*freq/(2*self.kT))/(1-math.exp(-self.h*self.c*100*freq/self.kT))
            self.vib_ps *= ps
    # calculate rotational partition sum of one molecule
    def cal_rot_ps(self):
        if self.data['mom_inert'][0] > 0:  # non-linear molecules
            self.rot_ps = pow(self.data['sym_num'], -1)*pow(self.kT/(self.h*self.c), 1.5)*pow(math.pi, 0.5)
            for I in self.data['mom_inert']:
                rot_const = self.h/(8*math.pi**2*self.c*I*self.moi_to_si)
                self.rot_ps *= pow(rot_const, -0.5)
        else:  # linear molecules
            moi = self.data['mom_inert'][1]*self.moi_to_si
            self.rot_ps = 8*pow(math.pi, 2)*moi*self.kT*pow(self.data['sym_num'], -1)*pow(self.h, -2)  # eq 10.19 in Cramer
    # this translational partition sum only works for one molecule, for 1 mol need to divide by NA! to
    # account for the multiplicity of space occupation
    def cal_trans_ps(self):
        self.trans_ps = self.V*pow(2*math.pi*self.data['mol_mass']*self.amu_to_si*self.kT, 1.5)/pow(self.h, 3)

    # calculate total partition sum of one molecule
    def cal_total_ps(self):
        self.cal_vib_ps()
        self.cal_rot_ps()
        self.cal_trans_ps()
        self.total_ps = self.vib_ps*self.rot_ps*self.trans_ps
        if partition == True:
            print("Partition function: ", self.total_ps)

    # calculate vibrational energy of 1mol molecule, the expression is
    # E = \sum_i {hw_i/2*(e^{beta hw_i}+1)/(e^{beta hw_i}-1)}
    def cal_vib_E(self):
        self.vib_E = 0
        for freq in self.data['vib_freq']:
            self.vib_E += (self.h*self.c*100*freq/2)*(math.exp(self.h*self.c*100*freq/self.kT)+1)/(math.exp(self.h*self.c*100*freq/self.kT)-1)
        self.vib_E *= self.NA*self.Jtokcal*self.kcaltoHa

    # calculate vibrational energy of 1mol molecule using quasi_RRHO approximation
    def cal_vib_E_quasi_RRHO(self):
        self.vib_E_quasi_RRHO = 0
        v0 = float(100.0)
        # v0 = float(sys.argv[2])        # cut off frequency, set as 100 cm^-1
        v0 = float(args.cutoff)        # overwrite cut-off frequency
        for freq in self.data['vib_freq']:
            # energy of pure vibrational mode using quantum vibrational energy
            E_vib = (self.h*self.c*100*freq/2)*(math.exp(self.h*self.c*100*freq/self.kT)+1)/(math.exp(self.h*self.c*100*freq/self.kT)-1)
            E_vib *= self.NA*self.Jtokcal*self.kcaltoHa
            # energy of pure rotational mode using classical rotational energy
            E_rot = 0.5*self.kT*self.NA*self.Jtokcal*self.kcaltoHa
            # weighting parameter
            w = 1/(1+math.pow(v0/freq, 4))
            self.vib_E_quasi_RRHO += w*E_vib + (1-w)*E_rot

    # calculate rotational energy of 1mol molecule, E = 1.5kT
    def cal_rot_E(self):
        self.rot_E = 1.5 * self.kT
        self.rot_E *= self.NA*self.Jtokcal*self.kcaltoHa

    # calculate translational energy of 1mol molecule, E = 1.5kT
    def cal_trans_E(self):
        self.trans_E = 1.5 * self.kT
        self.trans_E *= self.NA*self.Jtokcal*self.kcaltoHa

    # calculate total energy of 1mol molecule
    def cal_total_E(self):
        self.cal_vib_E()
        self.cal_vib_E_quasi_RRHO()
        self.cal_rot_E()
        self.cal_trans_E()
        self.get_scf_energy()
        RT = self.NA * self.kT * self.Jtokcal*self.kcaltoHa
        self.total_E = self.vib_E + self.rot_E + self.trans_E + RT  # added RT to convert from U to H using IGM
        self.total_E_quasi_RRHO = self.vib_E_quasi_RRHO + self.rot_E + self.trans_E + RT  # added RT to convert from U to H using IGM
        self.scf_plus_h_corr = self.scf_energy + self.total_E
        self.scf_plus_q_h_corr = self.scf_energy + self.total_E_quasi_RRHO
        print('=========================================')
        print('      Enthalpy Components (Ha)       ')
        print('=========================================')
        if verbose == True:
            print('             H_VIB: {0:8.8f}'.format(self.vib_E))
            print('       qRRHO H_VIB: {0:8.8f}'.format(self.vib_E_quasi_RRHO))
            print('             H_ROT: {0:8.8f}'.format(self.rot_E))
            print('           H_TRANS: {0:8.8f}'.format(self.trans_E))
            print('                RT: {0:8.8f}'.format(RT))
            print('           H_TOTAL: {0:8.8f}'.format(self.total_E))
            print('       qRRHO TOTAL: {0:8.8f}'.format(self.total_E_quasi_RRHO))
        print('        SCF+H_corr: {0:8.8f}'.format(self.scf_plus_h_corr))
        print('      SCF+q_H_corr: {0:8.8f}'.format(self.scf_plus_q_h_corr))

    # calculate vibrational free energy of 1mol molecule
    def cal_vib_FreeEnergy(self):
        self.vib_FreeEnergy = -self.kT*math.log(self.vib_ps)*self.NA*self.Jtokcal*self.kcaltoHa

    # calculate vibrational free energy of 1mol molecule using quasi_RRHO approximation
    def cal_vib_FreeEnergy_quasi_RRHO(self):
        self.vib_FreeEnergy_quasi_RRHO = self.vib_E_quasi_RRHO - self.T*self.vib_S_quasi_RRHO

    # calculate rotational free energy of 1mol molecule
    def cal_rot_FreeEnergy(self):
        self.rot_FreeEnergy = -self.kT*math.log(self.rot_ps)*self.NA*self.Jtokcal*self.kcaltoHa

    # calculate translational free energy of 1mol molecule, including the 1/NA! term
    def cal_trans_FreeEnergy(self):
        self.trans_FreeEnergy = -self.kT*(self.NA*math.log(self.trans_ps)-self.NA*math.log(self.NA)+self.NA)*self.Jtokcal*self.kcaltoHa

    # calculate total free energy of 1mol molecule
    def cal_total_FreeEnergy(self):
        self.cal_vib_S_quasi_RRHO()
        self.cal_vib_FreeEnergy()
        self.cal_vib_FreeEnergy_quasi_RRHO()
        self.cal_rot_FreeEnergy()
        self.cal_trans_FreeEnergy()
        self.total_FreeEnergy = self.vib_FreeEnergy + self.rot_FreeEnergy + self.trans_FreeEnergy
        self.total_FreeEnergy_quasi_RRHO = self.vib_FreeEnergy_quasi_RRHO + self.rot_FreeEnergy + self.trans_FreeEnergy
        self.scf_plus_g_corr = self.scf_energy + self.total_FreeEnergy
        self.scf_plus_q_g_corr = self.scf_energy + self.total_FreeEnergy_quasi_RRHO
        print('=========================================')
        print('     Free Energy Components (Ha)   ')
        print('=========================================')
        if verbose == True:
            print('               VIB: {0:8.8f}'.format(self.vib_FreeEnergy))
            print('         qRRHO VIB: {0:8.8f}'.format(self.vib_FreeEnergy_quasi_RRHO))
            print('               ROT: {0:8.8f}'.format(self.rot_FreeEnergy))
            print('             TRANS: {0:8.8f}'.format(self.trans_FreeEnergy))
            print('             TOTAL: {0:8.8f}'.format(self.total_FreeEnergy))
            print('       qRRHO TOTAL: {0:8.8f}'.format(self.total_FreeEnergy_quasi_RRHO))
        print('      SCF + G_corr: {0:8.8f}'.format(self.scf_plus_g_corr))
        print('SCF + qRRHO_G_corr: {0:8.8f}'.format(self.scf_plus_q_g_corr))

    # calculate vibrational entropy of 1mol molecule
    def cal_vib_S(self):
        self.vib_S = (self.vib_E - self.vib_FreeEnergy)/self.T

    # calculate vibrational entropy of 1mol molecule using quasi_RRHO approximation
    def cal_vib_S_quasi_RRHO(self):
        self.vib_S_quasi_RRHO = 0
        v0 = float(100.0)
        # v0 = float(sys.argv[2])        # cut off frequency, set as 100 cm^-1
        v0 = float(args.cutoff)        # overwrite cut-off frequency
        Bav = self.moi_to_si*(self.data['mom_inert'][0]+self.data['mom_inert'][1]+self.data['mom_inert'][2])/3    # average moment of inertia
        for freq in self.data['vib_freq']:
            omega = 2*math.pi*self.c*100*freq   # angular frequency
            mu = self.hbar/(2*omega)
            mu_prime = mu*Bav/(mu+Bav)
            w = 1/(1+math.pow(v0/freq, 4))
            # entropy of pure vibrational mode using quantum vibrational entropy
            S_vib = self.NA*self.k*(self.hbar*omega/(self.kT*(math.exp(self.hbar*omega/self.kT)-1))\
                -math.log(1-math.exp(-self.hbar*omega/self.kT)))*self.Jtokcal*self.kcaltoHa
            # entropy of pure rotational mode using classical rotational entropy
            S_rot = self.NA*self.k*(0.5+0.5*math.log(8*math.pi**3*mu_prime*self.kT/self.h**2))*self.Jtokcal*self.kcaltoHa
            self.vib_S_quasi_RRHO += w*S_vib + (1-w)*S_rot

    # calculate rotational entropy of 1mol molecule
    def cal_rot_S(self):
        self.rot_S = (self.rot_E - self.rot_FreeEnergy)/self.T

    # calculate translational entropy of 1mol molecule
    def cal_trans_S(self):
        self.trans_S = (self.trans_E - self.trans_FreeEnergy)/self.T

    # calculate total entropy of 1mol molecule
    def cal_total_S(self):
        self.cal_vib_S()
        self.cal_rot_S()
        self.cal_trans_S()
        self.total_S = (self.total_E - self.total_FreeEnergy)/self.T
        self.total_S_quasi_RRHO = (self.total_E_quasi_RRHO - self.total_FreeEnergy_quasi_RRHO)/self.T
        print('=========================================')
        print('      Entropy Components (Ha)      ')
        print('=========================================')
        if verbose == True:
            print('             S_VIB: {0:8.8f}'.format(self.vib_S))
            print('       qRRHO S_VIB: {0:8.8f}'.format(self.vib_S_quasi_RRHO))
            print('             S_ROT: {0:8.8f}'.format(self.rot_S))
            print('           S_TRANS: {0:8.8f}'.format(self.trans_S))
        print('       T * S_TOTAL: {0:8.8f}'.format(self.total_S*self.T))
        print(' T * qRRHO S_TOTAL: {0:8.8f}'.format(self.total_S_quasi_RRHO*self.T))

# def main():
if __name__ == "__main__":

    args = get_args()
    verbose = args.verbose
    partition = args.partition

    for filename in args.filenames:
        test = QC_freq_to_FreeEnergy(filename)
        print("=========================================")
        print("Extracting thermodynamic data from: \n", filename)
        print("=========================================")
        test.get_scf_energy()
        print("  Final SCF energy:", test.scf_energy)
        print(f"Using v0 = {args.cutoff} cm^-1 and T = {args.temp} K")
        test.get_useful_data()
        test.print_useful_data()
        test.cal_total_E()
        test.cal_total_ps()
        test.cal_total_FreeEnergy()
        test.cal_total_S()
        print("=========================================")
        print("For pasting into Excel:")
        print("E+Hcorr, TS, E+Gcorr")
        print(test.scf_plus_q_h_corr, test.total_S_quasi_RRHO*test.T, test.scf_plus_q_g_corr)

# main()
