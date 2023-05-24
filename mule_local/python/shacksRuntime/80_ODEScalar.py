from mule.JobCompileOptions import *

class ODEScalar:

    def __init__(self):

        ## ODE Scalar parameters
        self.u0 = None
        self.param_a = None
        self.param_b = None

    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''
        if not 'runtime.ode_scalar' in filter_list:
            if self.u0 != None:
                uniqueIDStr += '_ode_u0'+str(self.u0)
            if self.param_a != None:
                uniqueIDStr += '_ode_param_a'+str(self.param_a)
            if self.param_b != None:
                uniqueIDStr += '_ode_param_b'+str(self.param_b)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''

        if self.u0 != None:
            retRuntimeOptionsStr += " --u0="+str(self.u0)
        if self.param_a != None:
            retRuntimeOptionsStr += " --param-a="+str(self.param_a)
        if self.param_b != None:
            retRuntimeOptionsStr += " --param-b="+str(self.param_b)

        return retRuntimeOptionsStr
