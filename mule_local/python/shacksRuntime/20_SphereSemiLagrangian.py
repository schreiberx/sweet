
from mule.JobCompileOptions import *

class SphereSemiLagrangian:

    def __init__(self):
        self.semi_lagrangian_approximate_sphere_geometry = None



    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''

        if self.semi_lagrangian_approximate_sphere_geometry != None:
            uniqueIDStr += '_spap'+str(self.semi_lagrangian_approximate_sphere_geometry)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''

        if self.semi_lagrangian_approximate_sphere_geometry != None:
            retRuntimeOptionsStr += ' --semi-lagrangian-approximate-sphere-geometry='+str(self.semi_lagrangian_approximate_sphere_geometry)

        return retRuntimeOptionsStr

    