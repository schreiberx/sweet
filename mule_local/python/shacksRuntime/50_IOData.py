from mule.JobCompileOptions import *


class IOData:
    def __init__(self):

        self.output_timestep_size = None
        self.output_filename = None
        self.output_file_mode = None
        self.output_time_scale = None
        self.output_time_scale_inv = None


    def load_from_dict(self, d):

        pass


    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        idstr = ''
        return idstr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''

        if self.output_timestep_size != None:
            retRuntimeOptionsStr += ' -o '+str(self.output_timestep_size)

        if self.output_filename != None:
            if self.output_filename == "":
                retRuntimeOptionsStr += ' --output-file-name=-'
            else:
                retRuntimeOptionsStr += ' --output-file-name='+self.output_filename

        if self.output_file_mode != None:
            retRuntimeOptionsStr += ' --output-file-mode='+self.output_file_mode

        if self.output_time_scale != None:
            retRuntimeOptionsStr += ' --output-time-scale='+str(self.output_time_scale)

        if self.output_time_scale_inv != None:
            retRuntimeOptionsStr += ' --output-time-scale-inv='+str(self.output_time_scale_inv)


        return retRuntimeOptionsStr
    
    
    