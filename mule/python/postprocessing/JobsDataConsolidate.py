#! /usr/bin/env python3

from mule.InfoError import *
from mule.postprocessing.JobData import *
from mule.postprocessing.JobsData import *
import copy
import glob
import os
import re
import sys
import shutil
import pickle
import mule.utils
import mule.postprocessing.utils



class JobsDataConsolidate(InfoError):
    """
    Consolidate data:

    This is required if e.g. error vs. timestep size should
    be plotted for different time integration methods
    """

    def __init__(
            self,
            jobs_data,
            verbosity : int = 0
    ):
        """
        Initialize consolidation with job data
        """
        self.jobs_data = jobs_data
        self.verbosity = verbosity


    def create_groups(
        self,
        group_identifiers : list,
    ):
        """
        Group together particular job data, e.g. jobs using the same time stepping method

        Parameters:
        -----------
            group_identifiers: list
                List with attributes of jobs which are grouped together
                This also assembles the column name
        """

        groups = {}

        for jobdir, job_data in self.jobs_data.get_flattened_data().items():
            group_attributes = []
            group_attributes_short = []

            # collect values of group attributes
            for i in group_identifiers:
                # Only consider group if it really exists
                if i in job_data:
                    group_attributes_short.append(str(job_data[i]))
                    group_attributes.append(i+'_'+str(job_data[i]))
                else:
                    print("Job data "+i+" doesn't exist. Ignoring!")

            # create string identifier
            group_attributes_id = "__".join(group_attributes)
            group_attributes_short_id = "__".join(group_attributes_short)

            # create new group if it doesn't exist
            if not group_attributes_short_id in groups:
                groups[group_attributes_short_id] = {}
                print(group_attributes_short_id)

            # append job data to group
            groups[group_attributes_short_id][jobdir] = job_data

        return groups


def JobsData_GroupsCleanupPostprocessed(
        job_groups,         # from JobsDataConsolidate.create_groups()
        tag_cleanup_info,
        pickle_file_default_prefix,
        pint = False,
        solution_type = "physical"
):
    """
    tag_cleanup_info: List of dictionary
        ref_file_startswith:    Only if the reference output data file starts with this string
        tag_src:    Tag to get value from
        tag_dst:    Set this tag to this value
    """

    # Iterate over all groups
    for key, group in job_groups.items():

        # Iterate over each job in the group
        for job_id, job_data in group.items():

            # Get list of reference files
            ref_output_files = mule.postprocessing.utils.get_job_output_files(job_data)

            # Iterate over all reference files
            for ref_output_file in ref_output_files:

                if solution_type == "physical" and "_spec" in ref_output_file and ".csv" in ref_output_file:
                    continue

                for ci in tag_cleanup_info:
                    ref_file_starts_with = ci['ref_file_starts_with']
                    tag_src = ci['tag_src']
                    tag_dst = ci['tag_dst']

                    ## conflicts e.g. between *prog_h* and *prog_h_pert* ???
                    ###if ref_output_file.startswith(ref_file_starts_with):
                    if ref_output_file.startswith(ref_file_starts_with + "_t"):

                        # Determine basename of reference file
                        ref_basename = mule.utils.remove_file_ending(ref_output_file)

                        # Basename of pickle file
                        pickle_basename = f"{pickle_file_default_prefix}{ref_basename}"

                        # index is created with [pickle basename of file].[tag in pickle file]
                        jindex = f"{pickle_basename}.{tag_src}"

                        value = job_data[jindex]

                        if pint:
                            idx = jindex.find("_iter")
                            niter = jindex[idx:idx + 8]
                            tag_dst += niter

                        job_data[tag_dst] = value

                        print(f"{tag_dst}: {value}")




class JobsData_GroupsPlottingScattered:
    def __init__(
            self,
            groups,
            x_attribute_name : str,
            y_attribute_name : str,
            placeholder = None,
            data_filter = None,
            sort_data = True,
            verbosity = 0,
            meta_attribute_name = None
    ):
        """
        Create a data structure which is suitable for plotting

        Parameters:
        -----------
            groups: dict
                Dictionary with groups

            x_attribute_key: str
                Row Key: attribute which should serve as primary key
                (e.g. number of cores, time step size)

            y_attribute_name: str
                Column data: attribute which should serve as data field

            placeholder: ...
                Placeholder for data which doesn't exist

            data_filter: lambda(x, y, jobdata)
                Return false if data should be filtered out

            sort_data: boolean
                Sort data according to x-values
                This is important for line-connected points

        Return:
        -------
            dict { job_dir : { x_values: list, y_values: list }}
        """

        self.verbosity = verbosity
        self.groups = groups

        plot_data = {}

        for group_id in sorted(self.groups):
            group_jobs = self.groups[group_id]

            x_values = []
            y_values = []
            if meta_attribute_name != None:
                meta_values = []

            for jobdir, jobdata in group_jobs.items():
                if not x_attribute_name in jobdata:
                    print("")
                    print("WARNING: No data for attribute "+str(x_attribute_name)+" found")
                    print("WARNING: Job: "+str(jobdir))
                    continue

                x = jobdata[x_attribute_name]

                if not y_attribute_name in jobdata:
                    # Provide the option to filter out 'None'
                    if data_filter != None:
                        if data_filter(x, placeholder, jobdata):
                            continue
                    y = placeholder
                else:
                    y = jobdata[y_attribute_name]
                    if data_filter != None:
                        if data_filter(x, y, jobdata):
                            continue

                meta_value = None
                if meta_attribute_name != None:
                    if meta_attribute_name in jobdata:
                        meta_value = jobdata[meta_attribute_name]

                x_values.append(x)
                y_values.append(y)

                if meta_attribute_name != None:
                    meta_values.append(meta_value)

            plot_data[group_id] = {
                'x_values': x_values,
                'y_values': y_values
            }

            if meta_attribute_name != None:
                plot_data[group_id]['meta_values'] = meta_values


        if sort_data:
            for group_id, data in plot_data.items():

                sort_key = [float(i) for i in data['x_values']]

                x = data['x_values']
                y = data['y_values']

                y = [y for _,y in sorted(zip(sort_key,y))]
                x = [x for _,x in sorted(zip(sort_key,x))]

                data['x_values'] = x
                data['y_values'] = y

        self.data = plot_data



    def print(self):
        for group_key in sorted(self.data):
            group_data = self.data[group_key]

            print("Group '"+group_key+"'")

            if 'label' in group_data:
                print(" label: '"+group_data['label']+"'")

            for x, y in zip(group_data['x_values'], group_data['y_values']):
                print("    "+str(x)+" -> "+str(y))
            print("")



    def write(self, filename):
        with open(filename, 'w') as f:
            for group_key in sorted(self.data):
                group_data = self.data[group_key]

                f.write("group\t"+group_key+"\n")
                if 'label' in group_data:
                    f.write("label\t"+group_data['label']+"\n")

                f.write("x_values\ty_values\n")
                for x, y in zip(group_data['x_values'], group_data['y_values']):
                    f.write(str(x)+"\t"+str(y)+"\n")
                f.write("\n")



    def get_data(
        self
    ):
        return self.data



    def get_data_float(
        self
    ):
        data_plotting = copy.deepcopy(self.data)
        for key, values in data_plotting.items():
            x = []
            y = []
            for (i, j) in zip(values['x_values'], values['y_values']):
                if j != None:
                    x.append(float(i))
                    y.append(float(j))

            values['x_values'] = x
            values['y_values'] = y

        return data_plotting
        



class JobsData_GroupsDataTable:
    def __init__(
        self,
        groups : dict,
        primary_key_attribute_name : str,
        data_attribute_name : str,
        placeholder = None,
        sort_data = True,
        data_filter = None,
        verbosity = 0
    ):

        """
        Create a table-like data structure

        Parameters:
        -----------
            groups: dict
                Dictionary with groups

            primary_attribute_key: str
                Row Key: attribute which should serve as primary key

            data_attribute_name : str:
                Column data: attribute which should serve as data field


        Return:
        -------
            Table with filled in group data


        Example:
        --------
            Group was generated with ['timestepping_method']

            primary_attribute_key: number_of_cores
            data_attribute_name: simulation_wallclock_time

            Result:
            Table which looks like this

            -    | val(attr1) | val(attr2) | ...
            ---------------------------------------
            1    |        1.1 |        2.1 | ...
            2    |        2.4 |        nan | ...
            7    |       None |        nan | ...
            10    |       ....
            40    |
            ^
            |
            primary key
        """

        self.verbosity = verbosity
        self.groups = groups

        #
        # Determine full set of primary keys in case that primary key is missing somewhere
        #
        row_keys = []
        for group_id, group_jobs in self.groups.items():
            for jobdir, jobdata in group_jobs.items():
                if not primary_key_attribute_name in jobdata:
                    print("")
                    print("WARNING: No data for attribute "+primary_key_attribute_name+" found")
                    print("WARNING: Job: "+jobdir)
                    continue

                primary_key = jobdata[primary_key_attribute_name]
                if not primary_key in row_keys:
                    row_keys.append(primary_key)

        # Convert to floating point
        # Here, we assume that it's numeric data
        #row_keys = [float(i) for i in row_keys]

        # Sort the primary keys
        if sort_data:
            row_keys.sort()

        # get column names
        col_keys = sorted(list(self.groups.keys()))

        # get dimensions of table
        ncols = len(self.groups)
        nrows = len(row_keys)

        # Create table data
        data = [[None for i in range(ncols+1)] for j in range(nrows+1)]

        for group_id, group_jobs in self.groups.items():
            col_key = group_id
            col_id = col_keys.index(col_key)

            for jobkey, jobdata in group_jobs.items():
                if self.verbosity > 5:
                    print("Job: "+jobkey)

                if not primary_key_attribute_name in jobdata:
                    print("")
                    print("WARNING: No data for attribute "+primary_key_attribute_name+" found")
                    print("WARNING: Job key: "+jobkey)
                    continue

                row_key = jobdata[primary_key_attribute_name]
                row_id = row_keys.index(row_key)

                if data[row_id+1][col_id+1] != None:
                    self.print_data_table(data)
                    print("")
                    print("ERROR: Duplicate entry detected")
                    print("ERROR: This typically happens if either")
                    print("ERROR:  a) Groups are colliding")
                    print("ERROR:  b) axis variables are incorrect")
                    print("")
                    raise Exception("Duplicate entry!")

                if not data_attribute_name in jobdata:
                    print("WARNING: attribute "+data_attribute_name+" not found")
                    print("Job key: "+jobkey)
                    #for key, value in jobdata.items():
                    #    print(" + "+key+": "+str(value))

                    # Ignore missing data, will be filled in by placeholder :-)
                    #raise Exception("attribute '"+data_attribute_name+"' not found")

                else:
                    x = row_key
                    y = jobdata[data_attribute_name]
                    if data_filter != None:
                        if data_filter(x, y, jobdata):
                            continue

                    data[row_id+1][col_id+1] = y

        data[0][0] = '-'
        #data[0][0] = primary_key_attribute_name+'\\'+data_attribute_name

        # Setup row labels (primary keys)
        for i in range(len(row_keys)):
            data[i+1][0] = row_keys[i]

        # Setup col labels
        for j in range(len(col_keys)):
            data[0][j+1] = col_keys[j]

        if placeholder != None:
            # Replace None fields with placeholder
            for j in range(1,len(row_keys)):
                for i in range(1,len(col_keys)):
                    if data[j][i] == None:
                        data[j][i] = placeholder

        self.data = data



    def create_data_table_float(
        self,
        group_identifiers : list,
        primary_key_attribute_name : str,
        data_attribute_name : str,
        placeholder = None,
        data_filter = None
    ):
        data = self.create_data_table(group_identifiers, primary_key_attribute_name, data_attribute_name, data_filter=data_filter)

        # Replace None fields with placeholder
        for j in range(1,len(data)):
            for i in range(1,len(data[0])):
                if data[j][i] == None:
                    data[j][i] = placeholder
                else:
                    data[j][i] = float(data[j][i])

        return data




    def print(self):
        for row in self.data:
            print("\t".join([str(d) for d in row]))



    def write(self, filename):
        with open(filename, 'w') as f:
            for row in self.data:
                f.write("\t".join([str(d) for d in row])+"\n")





class JobsData_DataTable:
    def __init__(
        self,
        jobs_data : dict,
        primary_key_attribute_name : str,
        data_attribute_list : list,
        placeholder = None,
        sort_data = True,
        data_filter = None,
        verbosity = 0,
    ):

        """
        Create a table-like data structure

        Parameters:
        -----------
            jobs_data:
                Dictionary with Job Data

            primary_attribute_key: str
                Row Key: attribute which should serve as primary key

            data_attribute_list : list
                Column data: list of attributes which should serve as data field

        Return:
        -------
            Table with information

        Example:
        --------
            primary_attribute_key = 'number_of_ranks'
            data_attribute_list = ['timestepping', 'rexi_timestepping']

            Table which looks like this
            -    |      foo |      bar | ...
            ---------------------------------------
            1    |      1.1 |      2.1 | ...
            2    |      2.4 |      nan | ...
            7    |     None |      nan | ...
            10    |     ....
            40    |
            ^
            |
            primary key
        """

        self.verbosity = verbosity
        jobs_data_flattened = jobs_data.get_flattened_data()

        #
        # Determine full set of primary keys in case that primary key is missing somewhere
        #
        row_keys = []
        for job_id, jobdata in jobs_data_flattened.items():
            if not primary_key_attribute_name in jobdata:
                print("")
                print("WARNING: No data for attribute "+primary_key_attribute_name+" found")
                print("WARNING: Job: "+jobdir)
                continue

            primary_key = jobdata[primary_key_attribute_name]
            if not primary_key in row_keys:
                row_keys.append(primary_key)

        # Sort the primary keys
        # Not sure if this is always a good idea, but it makes sense for plots with numerical values
        if sort_data:
            row_keys.sort()

        # get column names
        col_keys = data_attribute_list

        # get dimensions of table
        ncols = len(col_keys)
        nrows = len(row_keys)

        # Create table data
        data = [[None for i in range(ncols+1)] for j in range(nrows+1)]

        for jobkey, jobdata in jobs_data_flattened.items():
            if self.verbosity > 5:
                print("Job: "+jobkey)

            for col_id in range(len(col_keys)):

                if not primary_key_attribute_name in jobdata:
                    print("")
                    print("WARNING: No data for attribute "+primary_key_attribute_name+" found")
                    print("WARNING: Job key: "+jobkey)
                    continue

                row_key = jobdata[primary_key_attribute_name]
                row_id = row_keys.index(row_key)

                col_key = col_keys[col_id]

                if data[row_id+1][col_id+1] != None:
                    self.print_data_table(data)
                    print("")
                    print("ERROR: Duplicate entry detected")
                    print("ERROR: This typically happens if either")
                    print("ERROR:  a) Groups are colliding")
                    print("ERROR:  b) axis variables are incorrect")
                    print("")
                    raise Exception("Duplicate entry!")

                if not col_key in jobdata:
                    print("WARNING: attribute "+str(col_key)+" not found")
                    print("Job key: "+jobkey)
                    #for key, value in jobdata.items():
                    #    print(" + "+key+": "+str(value))

                    # Ignore missing data, will be filled in by placeholder :-)
                    #raise Exception("attribute '"+data_attribute_name+"' not found")

                else:
                    x = row_key
                    y = jobdata[col_key]
                    if data_filter != None:
                        if data_filter(x, y, jobdata):
                            continue

                    data[row_id+1][col_id+1] = y

        data[0][0] = '-'
        #data[0][0] = primary_key_attribute_name+'\\'+data_attribute_name

        # Setup row labels (primary keys)
        for i in range(len(row_keys)):
            data[i+1][0] = row_keys[i]

        # Setup col labels
        for j in range(len(col_keys)):
            data[0][j+1] = col_keys[j]

        if placeholder != None:
            # Replace None fields with placeholder
            for j in range(1,len(row_keys)):
                for i in range(1,len(col_keys)):
                    if data[j][i] == None:
                        data[j][i] = placeholder

        self.data = data


    def get_data(self):
        return self.data


    def get_data_float(
            self,
            placeholder = None
        ):

        data = copy.deepcopy(self.data)
        # Replace None fields with placeholder
        for j in range(1,len(data)):
            for i in range(1,len(data[0])):
                if self.data[j][i] == None:
                    if placeholder != None:
                        data[j][i] = placeholder
                else:
                    data[j][i] = float(data[j][i])

        return data



    def print(self):
        for row in self.data:
            print("\t".join([str(d) for d in row]))

    def write(self, filename):
        with open(filename, 'w') as f:
            for row in self.data:
                f.write("\t".join([str(d) for d in row])+"\n")



class JobsDataPInTConsolidate(InfoError):
    """
    Consolidate PInT data:

    Postprocessing tools allowing to plot PInT output data
    e.g. error vs iteration
    """

    def __init__(
            self,
            jobs_data,
            verbosity : int = 0
    ):
        """
        Initialize consolidation with job data
        """
        self.jobs_data = jobs_data
        self.verbosity = verbosity


    def create_copy_for_each_iteration(self):
        """
        For each job in jobs_data, create a copy for each iteration
        Include iteration as a field in jobgeneration.pickle
        """

        jobs = self.jobs_data.get_flattened_data()
        for job in jobs.keys():

            path = jobs[job]['jobgeneration.job_dirpath'];

            ## don't make copies of copies!
            if "_iter" in path:
                continue;

            ## find PInT output files in path
            list_files = glob.glob(path + "/*iter*.csv");

            ## find max iter
            max_iter = -1;
            for f in list_files:
                ff = os.path.basename(f).split("_iter")
                niter = int(ff[1].split(".csv")[0])

                max_iter = max(max_iter, niter)

            print(job, max_iter)

            ## modify jobegenration.pickle
            with open(path + "/jobgeneration.pickle", 'rb') as ff:
                # Pickle the 'data' dictionary using the highest protocol available.
                pickle_data = pickle.load(ff)

            orig_job_dirpath = pickle_data['jobgeneration']['job_dirpath'];
            orig_job_unique_id = pickle_data['jobgeneration']['job_unique_id'];

            if max_iter >= 0:
                for niter in range(max_iter + 1):
                    suffix = "_iter" + str(niter);
                    shutil.copytree(path, path + suffix, dirs_exist_ok = True);
                    pickle_data['jobgeneration']['job_dirpath'] = orig_job_dirpath + suffix
                    pickle_data['jobgeneration']['job_unique_id'] = orig_job_unique_id + suffix
                    pickle_data['runtime']['iteration'] = niter
                    with open(path + suffix + "/jobgeneration.pickle", 'wb') as ff:
                        pickle.dump(pickle_data, ff, pickle.HIGHEST_PROTOCOL)

