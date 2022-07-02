#! /usr/bin/env python3


import numpy as np
import mule_local.postprocessing.shtnsfiledata as shtnsfiledata


class SphereDataSpectral:

    def __init__(self, filename = None, setup_physical=True):

        self.data_physical = None
        self.data_spectral = None

        if filename != None:
            self.read_file(filename, setup_physical=setup_physical)

        pass


    def read_file(self, filename, setup_physical=True):
        input_file = filename
        f = open(input_file, 'rb')
        content = f.read()
        f.close()


        self.file_info = {}

        acc = ""
        line_nr = 0
        fin_detected = False

        i = 0
        for i in range(0, len(content)):
            c = content[i]

            if c == ord('\n'):
                if line_nr == 0:
                    if acc == "SWEET":
                        pass
                    else:
                        raise Exception("Header not detected, stopping")

                else:
                    s = "DATA_TYPE"
                    if acc[0:len(s)] == s:
                        self.file_info['data_type'] = acc[len(s)+1:]

                    s = "MODES_N_MAX"
                    if acc[0:len(s)] == s:
                        self.file_info['modes_n_max'] = int(acc[len(s)+1:])

                    s = "MODES_M_MAX"
                    if acc[0:len(s)] == s:
                        self.file_info['modes_m_max'] = int(acc[len(s)+1:])

                    s = "NUM_ELEMENTS"
                    if acc[0:len(s)] == s:
                        self.file_info['num_elements'] = int(acc[len(s)+1:])

                    s = "GRID_TYPE"
                    if acc[0:len(s)] == s:
                        self.file_info['grid_type'] = acc[len(s)+1:]

                    s = "FIN"
                    if acc[0:len(s)] == s:
                        #print("FIN detected")
                        fin_detected = True
                        i += 1
                        break

                acc = ""
                line_nr += 1

            else:
                acc += chr(c)

        if not fin_detected:
            raise Exception("FIN not detected in file")

        if self.file_info['data_type'] != "SH_DATA":
            raise Exception("Wrong data type "+file_info['data_type']+" in binary file")

        if 0:
            print("*"*80)
            for key, value in self.file_info.items():
                print(str(key)+": "+str(value))
            print("*"*80)

        data = content[i:]

        self.data_spectral = np.frombuffer(data, dtype=np.complex128)

        if setup_physical:
            sfd = shtnsfiledata.shtnsfiledata()
            sfd.setup(self.file_info)

            self.data_physical = sfd.spec2phys(self.data_spectral)

            # Deprecated
            self.data = self.data_physical


