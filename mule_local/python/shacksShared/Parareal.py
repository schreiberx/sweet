
from mule.JobCompileOptions import *

class Parareal:

    def __init__(self):
        
        ## parareal parameters
        self.parareal_enabled = 0
        self.parareal_coarse_slices = None
        self.parareal_convergence_threshold = -1
        self.parareal_verbosity = 0
        self.parareal_max_simulation_time = None
        self.parareal_coarse_timestepping_method = None
        self.parareal_coarse_timestepping_order = 1
        self.parareal_coarse_timestepping_order2 = 1
        self.parareal_coarse_timestep_size = -1;
        self.parareal_load_ref_csv_files = 0;
        self.parareal_path_ref_csv_files = "";
        self.parareal_load_fine_csv_files = 0;
        self.parareal_path_fine_csv_files = "";
        self.parareal_store_iterations = 1;
        self.parareal_spatial_coarsening = None;
        self.parareal_max_iter = None;


    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''
        
        if not 'runtime.parareal' in filter_list:
            if self.parareal_enabled:
                uniqueIDStr += '_PARAREAL'
                if not 'runtime.parareal_coarse_slices' in filter_list:
                    uniqueIDStr += '_par_'+str(self.parareal_coarse_slices)
                if not 'runtime.parareal_coarse_timestepping_method' in filter_list:
                    uniqueIDStr += '_ptsm_'+str(self.parareal_coarse_timestepping_method)
                if not 'runtime.parareal_coarse_timestep_size' in filter_list:
                    uniqueIDStr += '_pDt_'+str(self.parareal_coarse_timestep_size)
                if not 'runtime.parareal_store_iterations' in filter_list:
                    uniqueIDStr += '_pStore_'+str(self.parareal_store_iterations)
                if not 'runtime.parareal_spatial_coarsening' in filter_list:
                    uniqueIDStr += '_pSpc_'+str(self.parareal_spatial_coarsening)
                if not 'runtime.parareal_max_iter' in filter_list:
                    uniqueIDStr += '_pMaxIter_'+str(self.parareal_max_iter)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''
        
        ## Parareal parameters
        if self.parareal_enabled:
            getRuntimeOptions += " --parareal-enable=1"
            getRuntimeOptions += " --parareal-coarse-slices="+str(self.parareal_coarse_slices)
            getRuntimeOptions += " --parareal-convergence-threshold="+str(self.parareal_convergence_threshold)
            getRuntimeOptions += " --parareal-verbosity="+str(self.parareal_verbosity)
            getRuntimeOptions += " --parareal-max-simulation-time="+str(self.parareal_max_simulation_time)
            getRuntimeOptions += " --parareal-coarse-timestepping-method="+str(self.parareal_coarse_timestepping_method)
            getRuntimeOptions += " --parareal-coarse-timestepping-order="+str(self.parareal_coarse_timestepping_order)
            getRuntimeOptions += " --parareal-coarse-timestepping-order2="+str(self.parareal_coarse_timestepping_order2)
            getRuntimeOptions += " --parareal-coarse-timestep-size="+str(self.parareal_coarse_timestep_size);
            getRuntimeOptions += " --parareal-load-ref-csv-files="+str(self.parareal_load_ref_csv_files);
            getRuntimeOptions += " --parareal-path-ref-csv-files="+str(self.parareal_path_ref_csv_files);
            getRuntimeOptions += " --parareal-load-fine-csv-files="+str(self.parareal_load_fine_csv_files);
            getRuntimeOptions += " --parareal-path-fine-csv-files="+str(self.parareal_path_fine_csv_files);
            getRuntimeOptions += " --parareal-store-iterations="+str(self.parareal_store_iterations);
            getRuntimeOptions += " --parareal-spatial-coarsening="+str(self.parareal_spatial_coarsening);

            if self.parareal_max_iter != None:
                getRuntimeOptions += " --parareal-max-iter="+str(self.parareal_max_iter);


        return retRuntimeOptionsStr

    