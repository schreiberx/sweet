
from mule.JobCompileOptions import *

class Parareal:

    def __init__(self):
        
        # Parareal
        self.parareal = 'none'
        self.parareal_scalar = 'disable'
        self.parareal_plane = 'disable'
        self.parareal_sphere = 'disable'
        self.parareal_plane_swe = 'disable'
        self.parareal_plane_burgers = 'disable'


    def getSConsParams(self):
        retval = ''
        # Parareal
        retval += ' --parareal='+self.parareal
        retval += ' --parareal-scalar='+self.parareal_scalar
        retval += ' --parareal-plane='+self.parareal_plane
        retval += ' --parareal-sphere='+self.parareal_sphere
        retval += ' --parareal-plane-swe='+self.parareal_plane_swe
        retval += ' --parareal-plane-burgers='+self.parareal_plane_burgers

        return retval


    def sconsAddOptions(self, scons):

        scons.AddOption(    '--parareal',
                dest='parareal',
                type='choice',
                choices=['none', 'serial','mpi'],
                default='none',
                help='Enable Parareal (none, serial, mpi) [default: %default]\nOnly works, if Parareal is supported by the simulation'
        )
        self.parareal = scons.GetOption('parareal')

        scons.AddOption(    '--parareal-scalar',
                dest='parareal_scalar',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal for scalar problems (enable, disable) [default: %default]'
        )
        self.parareal_scalar = scons.GetOption('parareal_scalar')

        scons.AddOption(    '--parareal-plane',
                dest='parareal_plane',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal on the plane (enable, disable) [default: %default]'
        )
        self.parareal_plane = scons.GetOption('parareal_plane')

        scons.AddOption(    '--parareal-sphere',
                dest='parareal_sphere',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal on the sphere (enable, disable) [default: %default]'
        )
        self.parareal_sphere = scons.GetOption('parareal_sphere')

        scons.AddOption(    '--parareal-plane-swe',
                dest='parareal_plane_swe',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal for SWE on the plane (enable, disable) [default: %default]'
        )
        self.parareal_plane_swe = scons.GetOption('parareal_plane_swe')
    
        scons.AddOption(    '--parareal-plane-burgers',
                dest='parareal_plane_burgers',
                type='choice',
                choices=['enable', 'disable'],
                default='0',
                help='Enable Parareal for Burgers on the plane (enable, disable) [default: %default]'
        )
        self.parareal_plane_burgers = scons.GetOption('parareal_plane_burgers')


    def sconsValidateOptions(self):
        
        if not self.parareal == 'none':
            if self.program == "parareal_ode":
                self.parareal_scalar = 'enable';
            elif self.program == 'swe_plane' or self.program == 'burgers':
                self.parareal_plane = 'enable';
                if self.program == 'swe_plane':
                    self.parareal_plane_swe = 'enable';
                elif self.program == 'burgers':
                    self.parareal_plane_burgers = 'enable';
            elif self.program == 'swe_sphere':
                self.parareal_sphere = 'enable';
                
                
                
    def sconsAddFlags(self, env):
        
        if self.parareal == 'none':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL=0'])
        elif self.parareal == 'serial':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL=1'])
        elif self.parareal == 'mpi':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL=2'])
        else:
            raise Exception("Invalid option '"+str(self.parareal)+"' for parareal method")
        
        
        if self.parareal_scalar == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL_SCALAR=1'])
        if self.parareal_plane == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL_PLANE=1'])
        if self.parareal_sphere == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL_SPHERE=1'])
        if self.parareal_plane_swe == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL_PLANE_SWE=1'])
        if self.parareal_plane_burgers == 'enable':
            env.Append(CXXFLAGS=['-DSWEET_PARAREAL_PLANE_BURGERS=1'])

