#! /usr/bin/env python3


import numpy as np


class SphereData:
    """
    Read sphere data from arbitrary sources:

        * with .sweet ending

        * with .csv ending

        * in physical/spectral space

    Warning: Not everything is supported, yet
    """

    def __init__(self, filename = None, setup_physical=False, setup_spectral=False):

        self.original_data_space = None

        self.data_physical = None
        self.data_spectral = None

        self.setup_physical = setup_physical
        self.setup_spectral = setup_spectral

        if filename.endswith(".sweet"):
            self._read_file_sweet(filename)

        elif filename.endswith(".csv"):
            raise Exception("TODO")

        else:
            raise Exception("Unknown file ending")

        if setup_physical:
            # Convert to physical space
            if self.data_physical == None:
                self._spectral_to_physical()

        if setup_spectral:
            # Convert to spectral space
            if self.data_spectral == None:
                self._physical_to_spectral()


    def _spectral_to_physical(self):
        from mule.postprocessing.SphereDataOperators import SphereDataOperators

        ops = SphereDataOperators(file_info=self.file_info)
        self.data_physical = ops.spec2phys(self.data_spectral)

        self.file_info['lats'] = ops.lats
        self.file_info['lons'] = ops.lons

    def _physical_to_spectral(self):
        from mule.postprocessing.SphereDataOperators import SphereDataOperators

        raise Exception("TODO")

        ops = SphereDataOperators(file_info=self.file_info)
        self.data_spectral = ops.phys2spec(self.data_physical)

        self.file_info['lats'] = ops.lats
        self.file_info['lons'] = ops.lons


    def _read_file_sweet(self, filename, setup_physical=True):
        f = open(filename, 'rb')
        content = f.read()
        f.close()

        self.file_info = {}

        acc = ""
        line_nr = 0
        fin_detected = False

        self.original_data_space = None

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
            print("ERROR: Only binary data supported")
            raise Exception("Wrong data type "+file_info['data_type']+" in binary file")

        if 0:
            print("*"*80)
            for key, value in self.file_info.items():
                print(str(key)+": "+str(value))
            print("*"*80)

        data = content[i:]

        self.data_spectral = np.frombuffer(data, dtype=np.cdouble)
