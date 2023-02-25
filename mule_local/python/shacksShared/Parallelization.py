
from mule.JobCompileOptions import *

class Parallelization:

    def __init__(self):
        self.num_threads_space = None


    def load_from_dict(self, d):
        if 'num_threads_space' in d:
            self.num_threads_space = int(d['num_threads_space'])


    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''

        if not 'parallelization' in filter_list:
            if self.num_threads_space != None:
                uniqueIDStr += '_nts'+str(self.num_threads_space)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''
        
        if self.num_threads_space != None:
            ret_runtime_options += ' --num-threads-space='+str(self.num_threads_space)

        return retRuntimeOptionsStr

    