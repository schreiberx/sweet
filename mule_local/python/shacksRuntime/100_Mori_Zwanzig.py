from mule.JobCompileOptions import *

class Mori_Zwanzig:

    def __init__(self):
        self.MZ_S_geostrophic_min =  None
        self.MZ_S_geostrophic_max =  None
        self.MZ_S_gravity_west_min = None
        self.MZ_S_gravity_west_max = None
        self.MZ_S_gravity_east_min = None
        self.MZ_S_gravity_east_max = None
        self.MZ_F_geostrophic_min =  None
        self.MZ_F_geostrophic_max =  None
        self.MZ_F_gravity_west_min = None
        self.MZ_F_gravity_west_max = None
        self.MZ_F_gravity_east_min = None
        self.MZ_F_gravity_east_max = None

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
                if self.MZ_S_geostrophic_min != None:
                    uniqueIDStr += '_MZ_S_min_geostr'+str(self.MZ_S_geostrophic_min)
                if self.MZ_S_geostrophic_max != None:
                    uniqueIDStr += '_MZ_S_max_geostr'+str(self.MZ_S_geostrophic_max)
                if self.MZ_S_gravity_west_min != None:
                    uniqueIDStr += '_MZ_S_min_gw'+str(self.MZ_S_gravity_west_min)
                if self.MZ_S_gravity_west_max != None:
                    uniqueIDStr += '_MZ_S_max_gw'+str(self.MZ_S_gravity_west_max)
                if self.MZ_S_gravity_east_min != None:
                    uniqueIDStr += '_MZ_S_min_ge'+str(self.MZ_S_gravity_east_min)
                if self.MZ_S_gravity_east_max != None:
                    uniqueIDStr += '_MZ_S_max_ge'+str(self.MZ_S_gravity_east_max)

                if self.MZ_F_geostrophic_min != None:
                    uniqueIDStr += '_MZ_F_min_geostr'+str(self.MZ_F_geostrophic_min)
                if self.MZ_F_geostrophic_max != None:
                    uniqueIDStr += '_MZ_F_max_geostr'+str(self.MZ_F_geostrophic_max)
                if self.MZ_F_gravity_west_min != None:
                    uniqueIDStr += '_MZ_F_min_gw'+str(self.MZ_F_gravity_west_min)
                if self.MZ_F_gravity_west_max != None:
                    uniqueIDStr += '_MZ_F_max_gw'+str(self.MZ_F_gravity_west_max)
                if self.MZ_F_gravity_east_min != None:
                    uniqueIDStr += '_MZ_F_min_ge'+str(self.MZ_F_gravity_east_min)
                if self.MZ_F_gravity_east_max != None:
                    uniqueIDStr += '_MZ_F_max_ge'+str(self.MZ_F_gravity_east_max)

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
        retRuntimeOptionsStr += " --MZ-S-geostr-min="             + str(self.MZ_S_geostrophic_min)
        retRuntimeOptionsStr += " --MZ-S-geostr-max="             + str(self.MZ_S_geostrophic_max)
        retRuntimeOptionsStr += " --MZ-S-gw-min="                 + str(self.MZ_S_gravity_west_min)
        retRuntimeOptionsStr += " --MZ-S-gw-max="                 + str(self.MZ_S_gravity_west_max)
        retRuntimeOptionsStr += " --MZ-S-ge-min="                 + str(self.MZ_S_gravity_east_min)
        retRuntimeOptionsStr += " --MZ-S-ge-max="                 + str(self.MZ_S_gravity_east_max)

        retRuntimeOptionsStr += " --MZ-F-geostr-min="             + str(self.MZ_F_geostrophic_min)
        retRuntimeOptionsStr += " --MZ-F-geostr-max="             + str(self.MZ_F_geostrophic_max)
        retRuntimeOptionsStr += " --MZ-F-gw-min="                 + str(self.MZ_F_gravity_west_min)
        retRuntimeOptionsStr += " --MZ-F-gw-max="                 + str(self.MZ_F_gravity_west_max)
        retRuntimeOptionsStr += " --MZ-F-ge-min="                 + str(self.MZ_F_gravity_east_min)
        retRuntimeOptionsStr += " --MZ-F-ge-max="                 + str(self.MZ_F_gravity_east_max)

        retRuntimeOptionsStr += " --MZ-epsilon="                   + str(self.MZ_epsilon)
        retRuntimeOptionsStr += " --MZ-F="                         + str(self.MZ_F)
        retRuntimeOptionsStr += " --MZ-timestepping-method-P="     + str(self.MZ_timestepping_method_P)
        retRuntimeOptionsStr += " --MZ-timestepping-method-Q="     + str(self.MZ_timestepping_method_Q)
        retRuntimeOptionsStr += " --MZ-timestepping-order-P="      + str(self.MZ_timestepping_order_P)
        retRuntimeOptionsStr += " --MZ-timestepping-order-Q="      + str(self.MZ_timestepping_order_Q)
        retRuntimeOptionsStr += " --MZ-timestepping-order2-P="     + str(self.MZ_timestepping_order2_P)
        retRuntimeOptionsStr += " --MZ-timestepping-order2-Q="     + str(self.MZ_timestepping_order2_Q)


        return retRuntimeOptionsStr


