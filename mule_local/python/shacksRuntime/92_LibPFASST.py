
from mule.JobCompileOptions import *

class LibPFASST:

    def __init__(self):
        # Parameters for LibPFASST
        self.libpfasst_nlevels = None
        self.libpfasst_niters = None
        self.libpfasst_nnodes = None
        self.libpfasst_nsweeps = None
        self.libpfasst_nodes_type = None
        self.libpfasst_coarsening_multiplier = None
        self.libpfasst_use_rexi = None
        self.libpfasst_implicit_coriolis_force = None
        self.libpfasst_use_rk_stepper = None
        self.libpfasst_u2 = None
        self.libpfasst_u4 = None
        self.libpfasst_u6 = None
        self.libpfasst_u8 = None
        self.libpfasst_u_fields = None


    def load_from_dict(self, d):
        
        pass

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        idstr = ''
        
        # LibPFASST attributes
        if not 'runtime.libpfasst' in filter_list:
            if self.libpfasst_nlevels != None:
                idstr += '_pf_nlev'+str(self.libpfasst_nlevels)
            if self.libpfasst_niters != None:
                idstr += '_pf_nit'+str(self.libpfasst_niters)
            if self.libpfasst_nnodes != None:
                idstr += '_pf_nnod'+str(self.libpfasst_nnodes)
            if self.libpfasst_nsweeps != None:
                idstr += '_pf_nswps'+str(self.libpfasst_nsweeps)
            if self.libpfasst_nodes_type != None:
                idstr += '_'+str(self.libpfasst_nodes_type)
            if self.libpfasst_coarsening_multiplier != None:
                idstr += '_pf_coars'+str(self.libpfasst_coarsening_multiplier)
            if self.libpfasst_use_rexi != None:
                idstr += '_pf_rexi'+str(self.libpfasst_use_rexi)
            if self.libpfasst_implicit_coriolis_force != None:
                idstr += '_pf_icf'+str(self.libpfasst_implicit_coriolis_force)
            if self.libpfasst_use_rk_stepper != None:
                idstr += '_pf_rk'+str(self.libpfasst_use_rk_stepper)

        return idstr


    def getRuntimeOptions(self):
        retval = ''
        
        # TODO: think about a way to check that we're in LibPFASST mode s.th. we do not have to check every if statement
        if self.libpfasst_nlevels != None:
            retval += ' --libpfasst-nlevels='+str(self.libpfasst_nlevels)
        if self.libpfasst_niters != None:
            retval += ' --libpfasst-niters='+str(self.libpfasst_niters)
        if self.libpfasst_nnodes != None:
            retval += ' --libpfasst-nnodes='+str(self.libpfasst_nnodes)
        if self.libpfasst_nsweeps != None:
            retval += ' --libpfasst-nsweeps='+str(self.libpfasst_nsweeps)
        if self.libpfasst_nodes_type != None:
            retval += ' --libpfasst-nodes-type='+str(self.libpfasst_nodes_type)
        if self.libpfasst_coarsening_multiplier != None:
            retval += ' --libpfasst-coarsening-multiplier='+str(self.libpfasst_coarsening_multiplier)
        if self.libpfasst_use_rexi != None:
            retval += ' --libpfasst-use-rexi='+str(self.libpfasst_use_rexi)
        if self.libpfasst_implicit_coriolis_force != None:
            retval += ' --libpfasst-implicit-coriolis-force='+str(self.libpfasst_implicit_coriolis_force)
        if self.libpfasst_use_rk_stepper != None:
            retval += ' --libpfasst-use-rk-stepper='+str(self.libpfasst_use_rk_stepper)
        if self.libpfasst_u2 != None:
            retval += ' --libpfasst-u2='+str(self.libpfasst_u2)
        if self.libpfasst_u4 != None:
            retval += ' --libpfasst-u4='+str(self.libpfasst_u4)
        if self.libpfasst_u6 != None:
            retval += ' --libpfasst-u6='+str(self.libpfasst_u6)
        if self.libpfasst_u8 != None:
            retval += ' --libpfasst-u8='+str(self.libpfasst_u8)
        if self.libpfasst_u_fields != None:
            retval += ' --libpfasst-u-fields='+str(self.libpfasst_u_fields)


        return retval

    