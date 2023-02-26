
from mule.JobCompileOptions import *

class XBraid:

    def __init__(self):

        ## XBraid parameters
        self.xbraid_enabled = 0
        self.xbraid_max_levels = None
        self.xbraid_skip = None
        self.xbraid_min_coarse = None
        self.xbraid_nrelax = None
        self.xbraid_nrelax0 = None
        self.xbraid_tol = None
        self.xbraid_tnorm = None
        self.xbraid_cfactor = None
        self.xbraid_cfactor0 = None
        self.xbraid_max_iter = None
        self.xbraid_fmg = None
        self.xbraid_fmg_vcyc = None
        self.xbraid_res = None
        self.xbraid_storage = None
        self.xbraid_print_level = None
        self.xbraid_access_level = None
        self.xbraid_run_wrapper_tests = None
        self.xbraid_fullrnorm = None
        self.xbraid_use_seq_soln = None
        self.xbraid_use_rand = None
        self.xbraid_pt = None
        self.xbraid_timestepping_method = None
        self.xbraid_timestepping_order = None
        self.xbraid_timestepping_order2 = None
        self.xbraid_viscosity_order = "2"
        self.xbraid_viscosity_coefficient = "0"
        self.xbraid_verbosity = None
        self.xbraid_load_ref_csv_files = None
        self.xbraid_path_ref_csv_files = None
        self.xbraid_load_fine_csv_files = None
        self.xbraid_path_fine_csv_files = None
        self.xbraid_store_iterations = None
        self.xbraid_spatial_coarsening = None


    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''
        if not 'runtime.xbraid' in filter_list:
            if self.xbraid_enabled:
                uniqueIDStr += '_XBRAID'
                if self.xbraid_max_levels != None:
                    uniqueIDStr += '_xb_max_l'+str(self.xbraid_max_levels)
                if self.xbraid_skip!= None:
                    idstr += '_xb_skip'+str(self.xbraid_skip)
                if self.xbraid_min_coarse != None:
                    idstr += '_xb_min_c'+str(self.xbraid_min_coarse)
                if self.xbraid_nrelax != None:
                    idstr += '_xb_nrlx'+str(self.xbraid_nrelax)
                if self.xbraid_nrelax0 != None:
                    idstr += '_xb_nrlx0'+str(self.xbraid_nrelax0)
                if self.xbraid_tol != None:
                    idstr += '_xb_tol'+str(self.xbraid_tol)
                if self.xbraid_tnorm != None:
                    idstr += '_xb_tnorm'+str(self.xbraid_tnorm)
                if self.xbraid_cfactor != None:
                    idstr += '_xb_cfr'+str(self.xbraid_cfactor)
                if self.xbraid_cfactor0 != None:
                    idstr += '_xb_cfr0'+str(self.xbraid_cfactor0)
                if self.xbraid_max_iter != None:
                    idstr += '_xb_max_i'+str(self.xbraid_max_iter)
                if self.xbraid_fmg != None:
                    idstr += '_xb_fmg'+str(self.xbraid_fmg)
                if self.xbraid_fmg_vcyc != None:
                    idstr += '_xb_fmg_vcyc'+str(self.xbraid_fmg_vcyc)
                if self.xbraid_res != None:
                    idstr += '_xb_res'+str(self.xbraid_res)
                if self.xbraid_storage != None:
                    idstr += '_xb_stg'+str(self.xbraid_storage)
                if self.xbraid_print_level != None:
                    idstr += '_xb_prt_l'+str(self.xbraid_print_level)
                if self.xbraid_access_level != None:
                    idstr += '_xb_acc_l'+str(self.xbraid_access_level)
                if self.xbraid_run_wrapper_tests != None:
                    idstr += '_xb_wrap'+str(self.xbraid_run_wrapper_tests)
                if self.xbraid_fullrnorm != None:
                    idstr += '_xb_frnm'+str(self.xbraid_fullrnorm)
                if self.xbraid_use_seq_soln != None:
                    idstr += '_xb_seq'+str(self.xbraid_use_seq_soln)
                if self.xbraid_use_rand != None:
                    idstr += '_xb_rand'+str(self.xbraid_use_rand)
                if self.xbraid_pt != None:
                    idstr += '_xb_pt'+str(self.xbraid_pt)
                if self.xbraid_timestepping_method != None:
                    idstr += '_xb_tsm'+str(self.xbraid_timestepping_method)
                if self.xbraid_timestepping_order != None:
                    idstr += '_xb_tso'+str(self.xbraid_timestepping_order)
                if self.xbraid_timestepping_order != None:
                    idstr += '_xb_tso2'+str(self.xbraid_timestepping_order2)
                if self.xbraid_viscosity_order != None:
                    idstr += '_xb_viscorder'+str(self.xbraid_viscosity_order)
                if self.xbraid_viscosity_coefficient != None:
                    idstr += '_xb_visccoeff'+str(self.xbraid_viscosity_coefficient)
                if self.xbraid_verbosity != None:
                    idstr += '_xb_verb'+str(self.xbraid_verbosity)
                if self.xbraid_load_ref_csv_files != None:
                    idstr += '_xb_load_ref'+str(self.xbraid_load_ref_csv_files)
                if self.xbraid_path_ref_csv_files != None:
                    idstr += '_xb_path_ref'+str(self.xbraid_path_ref_csv_files)
                if self.xbraid_load_fine_csv_files != None:
                    idstr += '_xb_load_fine'+str(self.xbraid_load_fine_csv_files)
                if self.xbraid_path_fine_csv_files != None:
                    idstr += '_xb_path_fine'+str(self.xbraid_path_fine_csv_files)
                if self.xbraid_store_iterations != None:
                    idstr += '_xb_store_iterations'+str(self.xbraid_store_iterations)
                if self.xbraid_spatial_coarsening != None:
                    idstr += '_xb_spc'+str(self.xbraid_spatial_coarsening)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''
        
        ## XBraid parameters
        if self.xbraid_enabled:
            retRuntimeOptionsStr += " --xbraid-enable=1"
            retRuntimeOptionsStr += " --xbraid-max-levels="+str(self.xbraid_max_levels)
            retRuntimeOptionsStr += " --xbraid-skip="+str(self.xbraid_skip)
            retRuntimeOptionsStr += " --xbraid-min-coarse="+str(self.xbraid_min_coarse)
            retRuntimeOptionsStr += " --xbraid-nrelax="+str(self.xbraid_nrelax)
            retRuntimeOptionsStr += " --xbraid-nrelax0="+str(self.xbraid_nrelax0)
            retRuntimeOptionsStr += " --xbraid-tol="+str(self.xbraid_tol)
            retRuntimeOptionsStr += " --xbraid-tnorm="+str(self.xbraid_tnorm)
            retRuntimeOptionsStr += " --xbraid-cfactor="+str(self.xbraid_cfactor)
            retRuntimeOptionsStr += " --xbraid-cfactor0="+str(self.xbraid_cfactor0)
            retRuntimeOptionsStr += " --xbraid-max-iter="+str(self.xbraid_max_iter)
            retRuntimeOptionsStr += " --xbraid-fmg="+str(self.xbraid_fmg)
            retRuntimeOptionsStr += " --xbraid-fmg-vcyc="+str(self.xbraid_fmg_vcyc)
            retRuntimeOptionsStr += " --xbraid-res="+str(self.xbraid_res)
            retRuntimeOptionsStr += " --xbraid-storage="+str(self.xbraid_storage)
            retRuntimeOptionsStr += " --xbraid-print-level="+str(self.xbraid_print_level)
            retRuntimeOptionsStr += " --xbraid-access-level="+str(self.xbraid_access_level)
            retRuntimeOptionsStr += " --xbraid-run-wrapper-tests="+str(self.xbraid_run_wrapper_tests)
            retRuntimeOptionsStr += " --xbraid-fullrnorm="+str(self.xbraid_fullrnorm)
            retRuntimeOptionsStr += " --xbraid-use-seq-soln="+str(self.xbraid_use_seq_soln)
            retRuntimeOptionsStr += " --xbraid-use-rand="+str(self.xbraid_use_rand)
            retRuntimeOptionsStr += " --xbraid-pt="+str(self.xbraid_pt)
            retRuntimeOptionsStr += " --xbraid-timestepping-method="+str(self.xbraid_timestepping_method)
            retRuntimeOptionsStr += " --xbraid-timestepping-order="+str(self.xbraid_timestepping_order2)
            retRuntimeOptionsStr += " --xbraid-timestepping-order2="+str(self.xbraid_timestepping_order2)
            retRuntimeOptionsStr += " --xbraid-viscosity-order="+str(self.xbraid_viscosity_order)
            retRuntimeOptionsStr += " --xbraid-viscosity-coefficient="+str(self.xbraid_viscosity_coefficient)
            retRuntimeOptionsStr += " --xbraid-verbosity="+str(self.xbraid_verbosity)
            retRuntimeOptionsStr += " --xbraid-load-ref-csv-files="+str(self.xbraid_load_ref_csv_files)
            retRuntimeOptionsStr += " --xbraid-path-ref-csv-files="+str(self.xbraid_path_ref_csv_files)
            retRuntimeOptionsStr += " --xbraid-load-fine-csv-files="+str(self.xbraid_load_fine_csv_files)
            retRuntimeOptionsStr += " --xbraid-path-fine-csv-files="+str(self.xbraid_path_fine_csv_files)
            retRuntimeOptionsStr += " --xbraid-store-iterations="+str(self.xbraid_store_iterations)
            retRuntimeOptionsStr += " --xbraid-spatial-coarsening="+str(self.xbraid_spatial_coarsening)

        return retRuntimeOptionsStr

    