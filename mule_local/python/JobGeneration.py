#! /usr/bin/env python3

import os
import stat

import hashlib

from mule.JobCompileOptions import *
from mule.JobRuntimeOptions import *

from mule.InfoError import *
from mule.JobPlatforms import *
from mule.JobPlatformResources import *
from mule.JobParallelization import *

__all__ = ['JobGeneration']

class JobGeneration(InfoError):

    def __init__(self, platform_id_override = None, dummy_init = False):
        self.init_phase = True

        InfoError.__init__(self, "JobGeneration")

        self.sweetroot = os.environ.get('MULE_SOFTWARE_ROOT')
        if self.sweetroot == None:
            self.error("Job environment variables not loaded!")

        # Setup all options
        self.compile : JobCompileOptions = JobCompileOptions(dummy_init=dummy_init)
        self.runtime : JobRuntimeOptions = JobRuntimeOptions(dummy_init=dummy_init)
        self.parallelization : JobParallelization = JobParallelization(dummy_init=dummy_init)

        self.platforms = JobPlatforms(platform_id_override, dummy_init=dummy_init)

        if dummy_init:
            self.platform_resources = JobPlatformResources()

            self.platform_id = "DUMMY"

        else:
            self.platform_functions = self.platforms.functions
            self.platform_resources = self.platform_functions.get_platform_resources()

            # Setup all hardware print_information which also ensures consistency
            self.platform_resources.setup()
            self.platform_resources.print()

            self.platform_id = self.platform_functions.get_platform_id()


        # Userdefined data for script header
        self.user_script_header = ''

        # Userdefined data for script footer
        self.user_script_footer = ''

        # Directly executed before the executable
        self.user_script_exec_prefix = ''

        # Directly executed after the executable
        self.user_script_exec_suffix = ''

        #
        # REXI parallelization with threading?
        # True deactivates OpenMP space parallelization and uses OMP parallel for over the REXI terms
        #self.rexi_thread_par = False

        # Should the compile command be part of the job script
        self.compilecommand_in_jobscript = True


        # List with filter names to reduce unique ID
        self.unique_id_filter = []

        # Unique ID of reference job to compare data to
        # This is extremely handy for postprocessing
        # This variable must be 'None' if this is a reference job
        self.reference_job_unique_id = None

        # Is this job a reference job?
        self.reference_job = False

        #
        # Unique ID of this job script
        # This is *** automatically generated *** if the job run.sh file is written!
        #
        self.job_unique_id = None

        # Path to directory
        self.job_dirpath = None

        #
        # p_ denote private variables
        # This should be changed
        #

        # Filepath to jobscript
        self.p_job_scriptpath = None

        # Directory where to execute job
        self.p_job_dirpath = None

        # Filepath to write output of jobscript to
        self.p_job_stdout_filepath = None
        self.p_job_stderr_filepath = None

        self.p_job_picklepath = None

        # Accumulator for compile commands
        self.p_compilecommands_accum = []


        self.init_phase = False


    def __setattr__(self, name, value):

        if name != 'init_phase':
            if not self.init_phase:
                if not name in self.__dict__:
                    raise Exception("Attribute '"+name+"' does not exist!")

        self.__dict__[name] = value



    def print(self):
        self.print_info("compilecommand_in_jobscript: "+str(self.compilecommand_in_jobscript))
        self.print_info("p_job_scriptpath: "+str(self.p_job_scriptpath))
        self.print_info("p_job_dirpath: "+str(self.p_job_dirpath))
        self.print_info("p_job_stdout_filepath: "+str(self.p_job_stdout_filepath))
        self.print_info("p_job_stderr_filepath: "+str(self.p_job_stderr_filepath))
        self.print_info("p_job_picklepath: "+str(self.p_job_picklepath))



    def setup_parallelization(self, parallelization_dim_list, **kwargs):
        """
        Setup the parallelization configuration

        Parameters
        ----------
        parallelization_dim_list: list
            list with JobParallelizationDimOptions with parallelization
            print_information along each dimension
        """
        self.parallelization.setup(parallelization_dim_list, self.platform_resources, **kwargs)



    def get_program_exec(self):
        """
        Return program executable and options

        This is accessed by the platform specific scripts to get the executable
        """
        return self.sweetroot+"/build/"""+self.compile.getProgramName()+""" """+self.runtime.getRuntimeOptions()



    def get_jobscript_content(self):
        """
        Generate the job script based on configuration
        """

        # Run Dummy setup in case that no setup was done so far!
        self.parallelization.dummy_setup_if_no_setup(self.platform_resources)

        if self.p_job_scriptpath == None:
            self.error("self.p_job_filepath == None")

        if self.p_job_dirpath == None:
            self.error("self.p_job_dirpath == None")

        if self.p_job_stdout_filepath == None:
            self.error("self.p_job_stdout_filepath == None")

        if self.p_job_stderr_filepath == None:
            self.error("self.p_job_stderr_filepath == None")



        self.compile.makeOptionsConsistent()

        content = ""

        # Setup job script generation
        self.platform_functions.jobscript_setup(self)

        # Now, assemble pieces together

        # Script header
        content += self.platform_functions.jobscript_get_header(self)
        content += self.user_script_header
        content += "# %SCRIPT_HEADER%\n"

        content += """

# Loading Job environment variables for currently active platform
cd \""""+self.sweetroot+"""\"
source ./local_software/env_vars.sh \""""+self.platforms.platform_id+"""\" || exit 1

"""
        #
        # Override compilers in environment variables
        #
        override_list = ['CC', 'CXX', 'F90', 'MPICC', 'MPICXX', 'MPIF90']

        for i in override_list:
            if 'MULE_'+i in os.environ:
                #print("INFO: Overriding environment variable "+i+"="+os.environ['MULE_'+i])
                content += "export "+i+"="+os.environ['MULE_'+i]+"\n"

        # Compile in main Job directory if requested
        c = self.platform_functions.jobscript_get_compile_command(self)
        if self.compilecommand_in_jobscript:
            content += c
        else:
            if c not in self.p_compilecommands_accum:
                self.p_compilecommands_accum.append(c)


        # Change to execution directory
        content += """

# chdir to execution directory
cd \""""+self.p_job_dirpath+"""\"

"""

        # EXEC prefix
        content += self.platform_functions.jobscript_get_exec_prefix(self)
        content += self.user_script_exec_prefix
        content += "# %SCRIPT_EXEC_PREFIX%\n"


        # EXEC
        content += self.platform_functions.jobscript_get_exec_command(self)

        # EXEC suffix
        content += "# %SCRIPT_EXEC_SUFFIX%\n"
        content += self.user_script_exec_suffix
        content += self.platform_functions.jobscript_get_exec_suffix(self)

        content += self.user_script_footer
        content += self.platform_functions.jobscript_get_footer(self)
        content += "# %SCRIPT_FOOTER%\n"

        #
        # Add return and exit 0
        #
        # return in case that it's sourced into someone's bash
        # exit 0 is required in case that the last program (e.g. a copy) failed
        #
        content += "return 2>/dev/null; exit 0"

        return content


    def get_compilecommands_content(self):
        """
        Return the content of a scriptfile to compile all programs for which the JobGeneration script was used for
        """
        content = """#! /bin/bash

# Loading Job environment variables
cd \""""+self.sweetroot+"""\"
source ./local_software/env_vars.sh \""""+os.path.normpath(self.platforms.platform_dirpath+"/env_vars.sh")+"""\" || exit 1

"""
        content += "\n".join(self.p_compilecommands_accum)

        # Empty accumulator
        self.p_compilecommands_accum = []
        return content

    def validate(self):
        if self.parallelization.max_wallclock_seconds != None and self.platform_resources.max_wallclock_seconds != None:
            if self.platform_resources.max_wallclock_seconds < self.parallelization.max_wallclock_seconds:
                self.info("platform_resources.max_wallclock_seconds: "+str(self.platform_resources.max_wallclock_seconds))
                self.info("parallelization.max_wallclock_seconds: "+str(self.parallelization.max_wallclock_seconds))
                self.error("Max. wallcock time exceeds platform's limit")


    def write_compilecommands(self, compilecommands_filename = None):
        if compilecommands_filename == None:
            compilecommands_filename = "./compile_platform_"+self.platforms.platform_id+".sh"

        compilecommands_filepath = os.path.abspath(compilecommands_filename)

        # Create directory for job script if it doesn't yet exist
        compilecommands_dirpath = os.path.dirname(compilecommands_filepath)

        # Generate script file content
        content = self.get_compilecommands_content()

        self.info("Compile commands file '"+compilecommands_filepath+"'")

        # Generate directory of scriptfile is not exist
        if not os.path.exists(compilecommands_dirpath):
            os.makedirs(compilecommands_dirpath)

        # Write scriptfile
        script_file = open(compilecommands_filepath, 'w')
        script_file.write(content)
        script_file.close()

        # Set permissions
        st = os.stat(compilecommands_filepath)
        os.chmod(compilecommands_filepath, st.st_mode | stat.S_IEXEC)



    def p_write_jobscript(self):

        # Validate before writing job script
        self.validate()

        # Generate script file content
        content = self.get_jobscript_content()

        # Write scriptfile
        script_file = open(self.p_job_scriptpath, 'w')
        self.info("Writing job file '"+self.p_job_scriptpath+"'")
        script_file.write(content)
        script_file.close()

        # Make script executable
        st = os.stat(self.p_job_scriptpath)
        os.chmod(self.p_job_scriptpath, st.st_mode | stat.S_IEXEC)



    def get_jobscript_directory(self):
        """
        Return the directory where the job data (including script) should be generated into
        """

        # Generate default jobname
        # Follow MULE's naming convention for job scripts
        if self.reference_job:
            i_job_dirpath = 'job_benchref_'+self.getUniqueID()
        else:
            i_job_dirpath = 'job_bench_'+self.getUniqueID()

        return i_job_dirpath



    def limit_filelength_with_hash(self, job_dirpath):
        # Don't allow more than 143 characters to support encrypted home directories
        max_length = 143

        l = len(job_dirpath)

        if l > max_length:
            m = hashlib.md5()
            m.update(job_dirpath.encode('utf-8'))

            # get hash
            h = m.hexdigest()
            job_dirpath = job_dirpath[0:max_length-len(h)-1]
            job_dirpath = job_dirpath + "_" + h

        return job_dirpath



    def gen_jobscript_directory(self, i_job_dirpath = None):
        """
        Generate a job script directory with
            run.sh:        Job script
            jobgeneration.pickle:    Information about job generation

        Parameters:
        --------------
        i_job_dirpath: str
            Directory to generate job script in.
            There is one directory per job!
        """

        if i_job_dirpath == None:
            i_job_dirpath = self.get_jobscript_directory()

        # IMPORTANT: Limit the max length of the job directory. Truncated part will be replaced with hash sum
        i_job_dirpath = self.limit_filelength_with_hash(i_job_dirpath)

        # Job directory
        self.p_job_dirpath = os.path.abspath(i_job_dirpath)

        # propagate to runtime script
        self.runtime.p_job_dirpath = self.p_job_dirpath

        self.job_dirpath = os.path.relpath(i_job_dirpath)

        # path to jobscript file
        self.p_job_scriptpath = self.p_job_dirpath+'/run.sh'

        # path to pickle file
        self.p_job_picklepath = self.p_job_dirpath+'/jobgeneration.pickle'

        # Determine output files
        self.p_job_stdout_filepath = self.p_job_dirpath+"/output.out"
        self.p_job_stderr_filepath = self.p_job_dirpath+"/output.out"

        # Generate full unfiltered unique ID
        self.job_unique_id = self.getUniqueID([])

        # Generate directory of scriptfile if it does not exist
        if not os.path.exists(self.p_job_dirpath):
            os.makedirs(self.p_job_dirpath)

        self.p_write_jobscript()
        self.save_file(self.p_job_picklepath)



    def gen_script(self, dirname, scriptname):
        """
        DEPRECATED
        DEPRECATED
        DEPRECATED
        """
        self.info("WARNING: gen_script(...) is deprecated, use gen_jobscript_directory(...)")

        self.write_jobscript(dirname+'/'+scriptname)



    def write_jobscript(self, jobscript_filepath):
        """
        DEPRECATED
        DEPRECATED
        DEPRECATED
        """
        self.info("WARNING: write_jobscript(...) THIS FUNCTION IS DEPRECATED")

        # Job directory
        self.p_job_dirpath = os.path.abspath(os.path.dirname(jobscript_filepath))
        self.job_dirpath = self.p_job_dirpath

        # path to jobscript file
        if os.path.basename(jobscript_filepath) != 'run.sh':
            raise Exception("Sorry, not supported anymore. Job script must be named 'run.sh'")

        self.p_job_scriptpath = self.p_job_dirpath+'/run.sh'

        # path to pickle file
        self.p_job_picklepath = self.p_job_dirpath+'/jobgeneration.pickle'

        # Determine output files
        self.p_job_stdout_filepath = self.p_job_dirpath+"/output.out"
        self.p_job_stderr_filepath = self.p_job_dirpath+"/output.err"

        # Generate directory of scriptfile if it does not exist
        if not os.path.exists(self.p_job_dirpath):
            os.makedirs(self.p_job_dirpath)

        self.p_write_jobscript()


    def getUniqueID(self, unique_id_filter = None):
        """
        Return a unique ID including *all* string and number attributes of this class
        """
        self.parallelization.dummy_setup_if_no_setup(self.platform_resources)

        if unique_id_filter == None:
            unique_id_filter = self.unique_id_filter

        unique_id = ''

        # First compile
        if 'compile' not in unique_id_filter:
            s = self.compile.getUniqueID(unique_id_filter)

            if s != '':
                if unique_id != '':
                    unique_id += '_'
                unique_id += s

        # Then runtime
        if 'runtime' not in unique_id_filter:
            s = self.runtime.getUniqueID(self.compile, unique_id_filter)
            if s != '':
                if unique_id != '':
                    unique_id += '_'
                unique_id += s

        # At the end the parallelization
        if 'parallelization' not in unique_id_filter:
            s = self.parallelization.getUniqueID(unique_id_filter)
            if s != '':
                if unique_id != '':
                    unique_id += '_'
                unique_id += s

        return unique_id



    def __get_sub_attributes_dict(self, obj, maxdepth=3):
        attr_dict = {}

        if maxdepth == 0:
            return attr_dict

        for key, a in obj.items():
            if isinstance(a, (float, int, str)):
                """
                Directly add this to the output list if the value is of type floating point, integer or string
                """
                attr_dict[key] = a


            elif isinstance(a, list):
                """
                Further trace down lists up to 'maxdepth'
                """
                for i in range(len(a)):
                    nattr = self.__get_sub_attributes_obj(a[i], maxdepth-1)

                    for key2, value in nattr.items():
                        attr_dict[key+"."+str(i)+"."+str(key2)] = value


            elif isinstance(a, dict):
                """
                Further trace down dictionaries up to 'maxdepth'
                """
                for key, aval in a.items():

                    if isinstance(aval, (float, int, str)):
                        attr_dict[attr+"."+str(key)] = value
                    else:
                        nattr = self.__get_sub_attributes_dict(aval, maxdepth-1)

                        for key2, value in nattr.items():
                            attr_dict[attr+"."+str(key)+"."+str(key2)] = value

        return attr_dict

    def __get_sub_attributes_obj(self, obj, maxdepth=2):
        attr_dict = {}

        if maxdepth == 0:
            return attr_dict

        obj_string_attributes = [a for a in dir(obj) if not a.startswith('__') and not callable(getattr(obj,a))]
        for attr in obj_string_attributes:
            a = getattr(obj, attr)

            if isinstance(a, (float, int, str)):
                """
                Directly add this to the output list if the value is of type floating point, integer or string
                """
                attr_dict[attr] = getattr(obj, attr)


            elif isinstance(a, list):
                """
                Further trace down lists up to 'maxdepth'
                """
                for i in range(len(a)):
                    nattr = self.__get_sub_attributes_obj(a[i], maxdepth-1)

                    for key, value in nattr.items():
                        attr_dict[attr+"."+str(i)+"."+str(key)] = value


            elif isinstance(a, dict):
                """
                Further trace down dictionaries up to 'maxdepth'
                """
                for key, aval in a.items():

                    if isinstance(aval, (float, int, str)):
                        attr_dict[attr+"."+str(key)] = value

                    elif isinstance(aval, (dict)):
                        nattr = self.__get_sub_attributes_dict(aval, maxdepth-1)

                        for key2, value in nattr.items():
                            attr_dict[attr+"."+str(key)+"."+str(key2)] = value

        return attr_dict


    def get_attributes_dict(self):
        attr_dict = {}
        attr_dict['jobgeneration'] = self.__get_sub_attributes_obj(self)
        attr_dict['compile'] = self.__get_sub_attributes_obj(self.compile)
        attr_dict['runtime'] = self.__get_sub_attributes_obj(self.runtime)
        attr_dict['parallelization'] = self.__get_sub_attributes_obj(self.parallelization)
        #attr_dict['platforms_platform'] = self.__get_sub_attributes_obj(self.platforms.platform)
        attr_dict['platform_resources'] = self.__get_sub_attributes_obj(self.platform_resources)
        attr_dict['platform_id'] = self.platform_id

        return attr_dict


    def save_file(self, filename):
        import pickle

        attr_dict = self.get_attributes_dict()

        with open(filename, 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(attr_dict, f, pickle.HIGHEST_PROTOCOL)


    def load_attributes_dict(self, filename):
        import pickle

        with open(filename, 'rb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            return pickle.load(f)


    def load_file(self, filename):
        def set_string_attributes(obj, attr_dict):
            for key, attr in attr_dict.items():
                setattr(obj, key, attr)

        data = self.load_attributes_dict(filename)

        set_string_attributes(self.compile, data['compile'])
        set_string_attributes(self.runtime, data['runtime'])
        set_string_attributes(self.parallelization, data['parallelization'])
        #set_string_attributes(self.platforms.platform, data['platforms_platform'])
        set_string_attributes(self.platform_resources, data['platform_resources'])
        self.platform_id = data['platform_id']



if __name__ == "__main__":
    p = JobGeneration()

    #s = p.get_jobscript_content()
    #p.info(s)

    s = p.getUniqueID()
    p.info(s)


    print("simtime: "+str(p.runtime.max_simulation_time))
    p.runtime.max_simulation_time = 666
    print("simtime: "+str(p.runtime.max_simulation_time))

    print("save_file()")
    p.save_file("/tmp/test.pickle")

    print("simtime: "+str(p.runtime.max_simulation_time))
    p.runtime.max_simulation_time = 123
    print("simtime: "+str(p.runtime.max_simulation_time))

    print("load_file()")
    p.load_file("/tmp/test.pickle")
    print("simtime: "+str(p.runtime.max_simulation_time))

    p.info("FIN")
