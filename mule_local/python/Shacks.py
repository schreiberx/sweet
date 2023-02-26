import os, sys

def getShacksDict(i_directory: str):
    """
    First, we will setup a list of all shack modules in the directory 'shacksShared'
    """
    import pkgutil
    from importlib.machinery import SourceFileLoader
    
    pkgpath = os.path.dirname(__file__)+"/"+i_directory
    shackFiles = [shackFile for _, shackFile, _ in pkgutil.iter_modules([pkgpath])]

    shacks = {}
    for shackFile in shackFiles:
        if shackFile[0] == '_':
            continue
        
        fullpath = pkgpath+"/"+shackFile+".py"
        module = SourceFileLoader("",fullpath).load_module()
        
        # Strip leading numbers and "_"
        shackName = shackFile[:]
        while shackName[0] in "0123456789_":
            shackName = shackName[1:]
        
        shacks[shackFile] = module.__dict__[shackName]()

    return shacks