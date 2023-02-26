
from mule.JobCompileOptions import *

class Polvani:
    def __init__(self):

        self.polvani_rossby = None
        self.polvani_froude = None


    def load_from_dict(self, d):
        if 'polvani_rossby' in d:
            self.polvani_rossby = float(d['polvani_rossby'])
            
        if 'polvani_froude' in d:
            self.polvani_froude = float(d['polvani_froude'])


    def getUniqueID(self,
                compileOptions : JobCompileOptions,
                filter_list : list = []
    ):
        uniqueIDStr = ''
        
        if not 'runtime.polvani' in filter_list:
            if self.polvani_rossby != None:
                uniqueIDStr += '_PR'+str(self.polvani_rossby)

            if self.polvani_froude != None:
                uniqueIDStr += '_PF'+str(self.polvani_froude)

        return uniqueIDStr


    def getRuntimeOptions(self):
        retRuntimeOptionsStr = ''
        
        if self.polvani_rossby != None:
            retRuntimeOptionsStr += ' --polvani-rossby='+str(self.polvani_rossby)

        if self.polvani_froude != None:
            retRuntimeOptionsStr += ' --polvani-froude='+str(self.polvani_froude)
        
        return retRuntimeOptionsStr

    