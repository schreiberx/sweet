from mule.JobCompileOptions import *

class Mori_Zwanzig:

    def __init__(self):
        self.MZ_SP_geostrophic_min =  None
        self.MZ_SP_geostrophic_max =  None
        self.MZ_SP_gravity_west_min = None
        self.MZ_SP_gravity_west_max = None
        self.MZ_SP_gravity_east_min = None
        self.MZ_SP_gravity_east_max = None
        self.MZ_SQ_geostrophic_min =  None
        self.MZ_SQ_geostrophic_max =  None
        self.MZ_SQ_gravity_west_min = None
        self.MZ_SQ_gravity_west_max = None
        self.MZ_SQ_gravity_east_min = None
        self.MZ_SQ_gravity_east_max = None
        self.MZ_FQ_geostrophic_min =  None
        self.MZ_FQ_geostrophic_max =  None
        self.MZ_FQ_gravity_west_min = None
        self.MZ_FQ_gravity_west_max = None
        self.MZ_FQ_gravity_east_min = None
        self.MZ_FQ_gravity_east_max = None

        self.MZ_epsilon = None
        self.MZ_F = None

        self.MZ_timestepping_method_P = None
        self.MZ_timestepping_method_Q = None
        self.MZ_timestepping_order_P = None
        self.MZ_timestepping_order_Q = None
        self.MZ_timestepping_order2_P = None
        self.MZ_timestepping_order2_Q = None

    def load_from_dict(self, d):
        return

    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''
        if not 'runtime.mori_zwanzig' in filter_list:
            if self.xbraid_enabled:
                uniqueIDStr += '_MZ'
                if self.MZ_SP_geostrophic_min != None:
                    uniqueIDStr += '_MZ_SP_min_geostr'+str(self.MZ_SP_geostrophic_min)
                if self.MZ_SP_geostrophic_max != None:
                    uniqueIDStr += '_MZ_SP_max_geostr'+str(self.MZ_SP_geostrophic_max)
                if self.MZ_SP_gravity_west_min != None:
                    uniqueIDStr += '_MZ_SP_min_gw'+str(self.MZ_SP_gravity_west_min)
                if self.MZ_SP_gravity_west_max != None:
                    uniqueIDStr += '_MZ_SP_max_gw'+str(self.MZ_SP_gravity_west_max)
                if self.MZ_SP_gravity_east_min != None:
                    uniqueIDStr += '_MZ_SP_min_ge'+str(self.MZ_SP_gravity_east_min)
                if self.MZ_SP_gravity_east_max != None:
                    uniqueIDStr += '_MZ_SP_max_ge'+str(self.MZ_SP_gravity_east_max)

                if self.MZ_SQ_geostrophic_min != None:
                    uniqueIDStr += '_MZ_SQ_min_geostr'+str(self.MZ_SQ_geostrophic_min)
                if self.MZ_SQ_geostrophic_max != None:
                    uniqueIDStr += '_MZ_SQ_max_geostr'+str(self.MZ_SQ_geostrophic_max)
                if self.MZ_SQ_gravity_west_min != None:
                    uniqueIDStr += '_MZ_SQ_min_gw'+str(self.MZ_SQ_gravity_west_min)
                if self.MZ_SQ_gravity_west_max != None:
                    uniqueIDStr += '_MZ_SQ_max_gw'+str(self.MZ_SQ_gravity_west_max)
                if self.MZ_SQ_gravity_east_min != None:
                    uniqueIDStr += '_MZ_SQ_min_ge'+str(self.MZ_SQ_gravity_east_min)
                if self.MZ_SQ_gravity_east_max != None:
                    uniqueIDStr += '_MZ_SQ_max_ge'+str(self.MZ_SQ_gravity_east_max)

                if self.MZ_FQ_geostrophic_min != None:
                    uniqueIDStr += '_MZ_FQ_min_geostr'+str(self.MZ_FQ_geostrophic_min)
                if self.MZ_FQ_geostrophic_max != None:
                    uniqueIDStr += '_MZ_FQ_max_geostr'+str(self.MZ_FQ_geostrophic_max)
                if self.MZ_FQ_gravity_west_min != None:
                    uniqueIDStr += '_MZ_FQ_min_gw'+str(self.MZ_FQ_gravity_west_min)
                if self.MZ_FQ_gravity_west_max != None:
                    uniqueIDStr += '_MZ_FQ_max_gw'+str(self.MZ_FQ_gravity_west_max)
                if self.MZ_FQ_gravity_east_min != None:
                    uniqueIDStr += '_MZ_FQ_min_ge'+str(self.MZ_FQ_gravity_east_min)
                if self.MZ_FQ_gravity_east_max != None:
                    uniqueIDStr += '_MZ_FQ_max_ge'+str(self.MZ_FQ_gravity_east_max)

                if self.MZ_epsilon != None:
                    uniqueIDStr += '_MZ_epsilon'+str(self.MZ_epsilon)
                if self.MZ_F != None:
                    uniqueIDStr += '_MZ_F'+str(self.MZ_F)
                if self.MZ_timestepping_method_P != None:
                    uniqueIDStr += '_MZ_tsm_P'+str(self.MZ_timestepping_method_P)
                if self.MZ_timestepping_method_Q != None:
                    uniqueIDStr += '_MZ_tsm_Q'+str(self.MZ_timestepping_method_Q)
                if self.MZ_timestepping_order_P != None:
                    uniqueIDStr += '_MZ_tso_P'+str(self.MZ_timestepping_order_P)
                if self.MZ_timestepping_order_Q != None:
                    uniqueIDStr += '_MZ_tso_Q'+str(self.MZ_timestepping_order_Q)
                if self.MZ_timestepping_order2_P != None:
                    uniqueIDStr += '_MZ_tso2_P'+str(self.MZ_timestepping_order2_P)
                if self.MZ_timestepping_order2_Q != None:
                    uniqueIDStr += '_MZ_tso2_Q'+str(self.MZ_timestepping_order2_Q)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''

        ## XBraid parameters
        retRuntimeOptionsStr += " --MZ-SP-geostr-min="             + str(self.MZ_SP_geostrophic_min)
        retRuntimeOptionsStr += " --MZ-SP-geostr-max="             + str(self.MZ_SP_geostrophic_max)
        retRuntimeOptionsStr += " --MZ-SP-gw-min="                 + str(self.MZ_SP_gravity_west_min)
        retRuntimeOptionsStr += " --MZ-SP-gw-max="                 + str(self.MZ_SP_gravity_west_max)
        retRuntimeOptionsStr += " --MZ-SP-ge-min="                 + str(self.MZ_SP_gravity_east_min)
        retRuntimeOptionsStr += " --MZ-SP-ge-max="                 + str(self.MZ_SP_gravity_east_max)

        retRuntimeOptionsStr += " --MZ-SQ-geostr-min="             + str(self.MZ_SQ_geostrophic_min)
        retRuntimeOptionsStr += " --MZ-SQ-geostr-max="             + str(self.MZ_SQ_geostrophic_max)
        retRuntimeOptionsStr += " --MZ-SQ-gw-min="                 + str(self.MZ_SQ_gravity_west_min)
        retRuntimeOptionsStr += " --MZ-SQ-gw-max="                 + str(self.MZ_SQ_gravity_west_max)
        retRuntimeOptionsStr += " --MZ-SQ-ge-min="                 + str(self.MZ_SQ_gravity_east_min)
        retRuntimeOptionsStr += " --MZ-SQ-ge-max="                 + str(self.MZ_SQ_gravity_east_max)

        retRuntimeOptionsStr += " --MZ-FQ-geostr-min="             + str(self.MZ_FQ_geostrophic_min)
        retRuntimeOptionsStr += " --MZ-FQ-geostr-max="             + str(self.MZ_FQ_geostrophic_max)
        retRuntimeOptionsStr += " --MZ-FQ-gw-min="                 + str(self.MZ_FQ_gravity_west_min)
        retRuntimeOptionsStr += " --MZ-FQ-gw-max="                 + str(self.MZ_FQ_gravity_west_max)
        retRuntimeOptionsStr += " --MZ-FQ-ge-min="                 + str(self.MZ_FQ_gravity_east_min)
        retRuntimeOptionsStr += " --MZ-FQ-ge-max="                 + str(self.MZ_FQ_gravity_east_max)

        retRuntimeOptionsStr += " --MZ-epsilon="                   + str(self.MZ_epsilon)
        retRuntimeOptionsStr += " --MZ-F="                         + str(self.MZ_F)
        retRuntimeOptionsStr += " --MZ-timestepping-method-P="     + str(self.MZ_timestepping_method_P)
        retRuntimeOptionsStr += " --MZ-timestepping-method-Q="     + str(self.MZ_timestepping_method_Q)
        retRuntimeOptionsStr += " --MZ-timestepping-order-P="      + str(self.MZ_timestepping_order_P)
        retRuntimeOptionsStr += " --MZ-timestepping-order-Q="      + str(self.MZ_timestepping_order_Q)
        retRuntimeOptionsStr += " --MZ-timestepping-order2-P="     + str(self.MZ_timestepping_order2_P)
        retRuntimeOptionsStr += " --MZ-timestepping-order2-Q="     + str(self.MZ_timestepping_order2_Q)


        return retRuntimeOptionsStr


