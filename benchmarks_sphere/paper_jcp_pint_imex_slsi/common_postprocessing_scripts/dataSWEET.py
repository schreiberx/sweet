import numpy as np
import matplotlib.pyplot as plt
import sys
from mule.SWEETRuntimeParametersScenarios import *
from mule.postprocessing.JobsData import *
from mule.postprocessing.JobsDataConsolidate import *
from common import *

class dataSWEET:

    def __init__(self, output_type, geometry):
        self.output_type = output_type;
        self.geometry = geometry;

        if self.geometry == "plane" or self.geometry == "scalar":
            self.timescale = 1;
        elif self.geometry == "sphere":
            self.timescale = 1. / (60. * 60.);
        else:
            raise Exception("Wrong geometry: ", self.geometry);

    ## get parameter value corresponding to given job
    def getValFromVar(self, var, job):
        return getValFromVar(var, job, self.jobs_data);

    def loadDataFromFile(self, filename):

        if self.output_type == "csv":
            data = self.loadDataFromFileCSV(filename + ".csv");
        elif self.output_type == "bin":
            data = self.loadDataFromFileBin(filename + ".sweet");
        else:
            sys.exit("Unknown output_type: " + self.output_type);

        return data


    ## read CSV files
    def loadDataFromFileCSV(self, filename):

        if self.geometry == "scalar":
            try:
                data = {};
                lines = [line.rstrip() for line in open(filename)];
                for i in range(3, len(lines)):
                    spl = lines[i].split();
                    try:
                        data[i-3] = np.array(float(spl[0]));
                    except:
                        try:
                            idx = spl[0].find(",");
                            real = float(spl[0][1:idx]);
                            imag = float(spl[0][idx+1:-1]);
                            data[i-3] = np.array(complex(real, imag));
                        except:
                            data[i-3] = None;
            except:
                data = None;
            return data;

        try:
            data = np.loadtxt(filename)
            ## remove first row and columns (corresponding to lat/lon values)
            if self.geometry == "sphere":
                data = data[1:, 1:]
        except:
            data = None;

        return data

    ## read BIN files (return solution converted to physical space)
    def loadDataFromFileBin(self, filename):
        f = open(filename, 'rb')
        content = f.read()
        f.close()

        numbers = re.findall("\d+", filename)
        time = int(numbers[-2])/24

        file_info = {}

        acc = ""
        line_nr = 0
        fin_detected = False

        i = 0
        for i in range(0, len(content)):
            c = content[i]

            if c == ord('\n'):
                if line_nr == 0:
                    if acc == "SWEET":
                        print("SWEET header detected")
                    else:
                        raise Exception("Header not detected, stopping")

                else:
                    s = "DATA_TYPE"
                    if acc[0:len(s)] == s:
                        file_info['data_type'] = acc[len(s)+1:]

                    s = "MODES_N_MAX"
                    if acc[0:len(s)] == s:
                        file_info['modes_n_max'] = int(acc[len(s)+1:])

                    s = "MODES_M_MAX"
                    if acc[0:len(s)] == s:
                        file_info['modes_m_max'] = int(acc[len(s)+1:])

                    s = "NUM_ELEMENTS"
                    if acc[0:len(s)] == s:
                        file_info['num_elements'] = int(acc[len(s)+1:])

                    s = "GRID_TYPE"
                    if acc[0:len(s)] == s:
                        file_info['grid_type'] = acc[len(s)+1:]

                    s = "FIN"
                    if acc[0:len(s)] == s:
                        print("FIN detected")
                        fin_detected = True
                        i += 1
                        break

                acc = ""
                line_nr += 1

            else:
                acc += chr(c)

        if not fin_detected:
            raise Exception("FIN not detected in file")

        if file_info['data_type'] != "SH_DATA":
            raise Exception("Wrong data type "+file_info['data_type']+" in binary file")

        print("*"*80)
        for key, value in file_info.items():
            print(str(key)+": "+str(value))
        print("*"*80)

        data = content[i:]
        print("BINARY DATA LEN: "+str(len(data)))

        s = SHTNS_data()
        data_spec = s.setup(file_info, data)
        data_phys = s.spec2phys(data_spec)

        return data_phys;

    def findFineSimulation(self, pint_type):

        if not (pint_type == "parareal" or pint_type == "xbraid"):
            sys.exit("Wrong pint_type:", pint_type);

        ###########print("Finding fine simulation");
        ###########key_fine = "";
        ###########for job in self.jobs_data.keys():
        ###########    if not self.jobs_data[job]["runtime." + pint_type + "_enabled"]:
        ###########        if not key_fine == "":
        ###########            sys.exit("There are more than one fine simulation!")
        ###########        key_fine = job;
        ###########if key_fine == "":
        ###########    sys.exit("No fine simulation has been found!");

        print("Finding fine and ref simulation");
        keys = [];
        for job in self.jobs_data.keys():
            if not self.jobs_data[job]["runtime." + pint_type + "_enabled"]:
                keys.append(job);
        if len(keys) == 0:
            sys.exit("No fine or ref simulations have been found")
        elif len(keys) > 2:
            sys.exit("More than one fine and ref simulations have been found")
        elif len(keys) == 1:
            key_fine = keys[0]
            key_ref = keys[0]
        else:
            dt1 = self.jobs_data[keys[0]]["runtime.timestep_size"]
            dt2 = self.jobs_data[keys[1]]["runtime.timestep_size"]
            if (dt1 < dt2):
                key_ref = keys[0]
                key_fine = keys[1]
            else:
                key_ref = keys[1]
                key_fine = keys[0]

        return key_fine, key_ref;

    ## create file with jobs details (values of given list of parameters)
    def storeJobsDetails(self, list_params):

        f = open("job_details", "w");
        for job in self.jobs_data.keys():
            f.write("\n" + getJobPath(self.jobs_data, job) + "\n");
            for param in list_params:
                try:
                    val = self.jobs_data[job][param];
                    f.write(param + "\t\t\t\t\t" + str(val) + "\n");
                except:
                    pass;
            f.write("\n");
        f.close();


    ## read all jobs in folder, except if a list is given
    def readAllData(self, prefix_files, t_input = None):

        self.jobs_solution = {};

        ## read simulations
        self.jobs_data = JobsData('job_bench*', verbosity=0).get_flattened_data();

        for key in self.jobs_data.keys():

            path = self.jobs_data[key]["jobgeneration.p_job_dirpath"];
            path = os.path.basename(path);
            dt = self.jobs_data[key]["runtime.timestep_size"];
            dt_output = self.jobs_data[key]["runtime.output_timestep_size"];
            Tmax = self.jobs_data[key]["runtime.max_simulation_time"];
            parareal_enabled = self.jobs_data[key]["runtime.parareal_enabled"];
            xbraid_enabled = self.jobs_data[key]["runtime.xbraid_enabled"];

            parareal_enabled = self.jobs_data[key]["runtime.parareal_enabled"];
            xbraid_enabled = self.jobs_data[key]["runtime.xbraid_enabled"];
            self.jobs_solution[key] = {};

            ## read reference / fine solution
            if (not parareal_enabled) and (not xbraid_enabled):
                nb_timesteps = int(Tmax / dt);
                for it in range(0, nb_timesteps + 1):
                    if not (t_input is None):
                        if abs(t_input - it * Tmax / nb_timesteps ) > 1e-10:
                            continue;
                    tt_str = getFormattedDate(float(it * Tmax / nb_timesteps), self.timescale);
                    filename = path + "/" + prefix_files + tt_str;
                    self.jobs_solution[key][tt_str] = self.loadDataFromFile(filename);
            ## read pint solution
            else:
                if parareal_enabled:
                    parareal_coarse_slices = self.jobs_data[key]["runtime.parareal_coarse_slices"];
                    for it in range(0, parareal_coarse_slices + 1):
                        tt_str = self.getFormattedDate(float(it * Tmax / parareal_coarse_slices), self.timescale);
                        self.jobs_solution[key][tt_str] = {};
                        for niter in range(min(it + 1, parareal_coarse_slices + 1)):
                            if it == 0:
                                filename = path + "/" + prefix_files + tt_str;
                            else:
                                filename = path + "/" + prefix_files + tt_str + "_iter" + getFormattedIteration(niter);
                            self.jobs_solution[key][tt_str][niter] = self.loadDataFromFile(filename);

                if xbraid_enabled:
                   xbraid_max_iter = self.jobs_data[key]["runtime.xbraid_max_iter"];
                   nb_timesteps = int(Tmax / dt_output);
                   for it in range(0, nb_timesteps + 1):
                       if not (t_input is None):
                           if abs(t_input - it * dt_output ) > 1e-10:
                               continue;
                       tt_str = getFormattedDate(float(it * Tmax / nb_timesteps), self.timescale);
                       self.jobs_solution[key][tt_str] = {};
                       for niter in range(xbraid_max_iter):
                           if it == 0:
                               filename = path + "/" + prefix_files + tt_str + "_iter" + getFormattedIteration(niter);
                           else:
                               filename = path + "/" + prefix_files + tt_str + "_iter" + getFormattedIteration(niter);
                           self.jobs_solution[key][tt_str][niter] = self.loadDataFromFile(filename);


        return self.jobs_data, self.jobs_solution;


    ## some Parareal and Xbraid parameters are in the form x,y
    ## where x and y correspond respectively to the fine and the coarse levels
    ## This function splits the param "param" into "param_fine" and "param_coarse"
    ## If there are more than two values (MGRIT with nlevels > 2), then create "param_coarse2" etc.
    def splitPinTParameters(self, pint_type):

        list_vars = [
                        "runtime." + pint_type + "_timestepping_method",
                        "runtime." + pint_type + "_timestepping_order",
                        "runtime." + pint_type + "_timestepping_order2",
                        "runtime." + pint_type + "_viscosity_order",
                        "runtime." + pint_type + "_viscosity_coefficient"
                    ]

        for job in self.jobs_data.keys():

            ###for var in self.jobs_data[job].copy().keys():

            for var in list_vars:

                if not (var in self.jobs_data[job].copy().keys()):
                    continue

                if not pint_type in var:
                    continue
                val = self.jobs_data[job][var]
                ####if not "," in str(val):
                ####    continue
                if type(val) is str:
                    spl = val.split(",")
                else:
                    spl = [str(val)]
                try:
                    if 'viscosity_order' in var:
                        for i in range(len(spl)):
                            spl[i] = int(spl[i])
                    else:
                        for i in range(len(spl)):
                            spl[i] = float(spl[i])
                except:
                    pass

                self.jobs_data[job][var + "_fine"] = spl[0]
                if len(spl) > 1:
                    self.jobs_data[job][var + "_coarse"] = spl[1]
                    self.jobs_data[job][var + "_coarse1"] = spl[1]
                    for i in range(2, len(spl)):
                        self.jobs_data[job][var + "_coarse" + str(i)] = spl[i]
                else:
                    self.jobs_data[job][var + "_coarse"] = spl[0]

                ## repeat values for coarse level if a single value (= fine) was given
                if len(spl) == 1:
                    for i in range(1, self.jobs_data[job]["runtime.xbraid_max_levels"]):
                        self.jobs_data[job][var + "_coarse" + str(i)] = spl[0]

                ## repeat values for coarse level if a single coarse value was given
                if len(spl) == 2:
                    for i in range(2, self.jobs_data[job]["runtime.xbraid_max_levels"]):
                        self.jobs_data[job][var + "_coarse" + str(i)] = spl[1]

                ## add dummy values
                for i in range(self.jobs_data[job]["runtime.xbraid_max_levels"], len(spl) + 5):
                    self.jobs_data[job][var + "_coarse" + str(i)] = "--"

                ###self.jobs_data[job][var + "_allcoarse"] = spl[1]
        return self.jobs_data

    def getListJobs(self):
        return list(self.jobs_data.keys());
