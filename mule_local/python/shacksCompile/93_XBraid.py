
from mule.JobCompileOptions import *

class XBraid:

    def __init__(self):
        # XBraid
        self.xbraid = 'none'
        self.xbraid_scalar = 'disable'
        self.xbraid_plane = 'disable'
        self.xbraid_sphere = 'disable'
        self.xbraid_plane_swe = 'disable'
        self.xbraid_plane_burgers = 'disable'


    def getSConsParams(self):
        retval = ''

        # XBraid
        retval += ' --xbraid='+self.xbraid
        retval += ' --xbraid-scalar='+self.xbraid_scalar
        retval += ' --xbraid-plane='+self.xbraid_plane
        retval += ' --xbraid-sphere='+self.xbraid_sphere
        retval += ' --xbraid-plane-swe='+self.xbraid_plane_swe
        retval += ' --xbraid-plane-burgers='+self.xbraid_plane_burgers


        return retval


    def sconsAddOptions(self, scons):

        scons.AddOption(    '--xbraid',
                dest='xbraid',
                type='choice',
                choices=['none', 'mpi'],
                default='none',
                help='Enable XBBraid (none,  mpi) [default: %default]\nOnly works, if XBraid is supported by the simulation'
        )
        self.xbraid = scons.GetOption('xbraid')

        scons.AddOption(    '--xbraid-scalar',
                dest='xbraid_scalar',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for scalar problems (enable, disable) [default: %default]'
        )
        self.xbraid_scalar = scons.GetOption('xbraid_scalar')

        scons.AddOption(    '--xbraid-plane',
                dest='xbraid_plane',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid on the plane (enable, disable) [default: %default]'
        )
        self.xbraid_plane = scons.GetOption('xbraid_plane')

        scons.AddOption(    '--xbraid-sphere',
                dest='xbraid_sphere',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid on the sphere (enable, disable) [default: %default]'
        )
        self.xbraid_sphere = scons.GetOption('xbraid_sphere')

        scons.AddOption(    '--xbraid-plane-swe',
                dest='xbraid_plane_swe',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for SWE on the plane (enable, disable) [default: %default]'
        )
        self.xbraid_plane_swe = scons.GetOption('xbraid_plane_swe')

        scons.AddOption(    '--xbraid-plane-burgers',
                dest='xbraid_plane_burgers',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable XBraid for Burgers on the plane (enable, disable) [default: %default]'
        )
        self.xbraid_plane_burgers = scons.GetOption('xbraid_plane_burgers')


    def sconsValidateOptions(self):
        
        if not self.xbraid == 'none':
            if self.program == "parareal_ode":
                self.xbraid_scalar = 'enable';
            elif self.program == 'swe_plane' or self.program == 'burgers':
                self.xbraid_plane = 'enable';
                if self.program == 'swe_plane':
                    self.xbraid_plane_swe = 'enable';
                elif self.program == 'burgers':
                    self.xbraid_plane_burgers = 'enable';
            elif self.program == 'swe_sphere':
                self.xbraid_sphere = 'enable';
                
                
    def sconsAddFlags(self, env):
        if self.xbraid == 'none':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID=0'])
        elif self.xbraid == 'mpi':
            env.Append(CXXFLAGS=['-Ilocal_software/local/include/xbraid'])
            env.Append(LIBS=['braid'])
            env.Append(CXXFLAGS=['-DSWEET_XBRAID=1'])
        else:
            raise Exception("Invalid option '"+str(self.xbraid)+"' for XBraid")
        
        
        if self.xbraid_scalar == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_SCALAR=1'])
        if self.xbraid_plane == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_PLANE=1'])
        if self.xbraid_sphere == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_SPHERE=1'])
        if self.xbraid_plane_swe == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_PLANE_SWE=1'])
        if self.xbraid_plane_burgers == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_XBRAID_PLANE_BURGERS=1'])
        

