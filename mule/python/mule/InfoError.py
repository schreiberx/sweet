#! /usr/bin/env python3


__all__ = ['InfoError']

#
# Usage: e.g.
#
# class  pickle_PlaneDataPhysicalDiff(InfoError):
# def __init__(self):
#     InfoError.__init__(self, "pickle_PlaneDataPhysicalDiff")
#     ...
#


class InfoError:
    def __init__(self, i_prefix):
        InfoError.setup(self, i_prefix)

    def setup(self, i_prefix):
        self.prefix = i_prefix

    def info(self, i_str):
        o = "INFO ["+self.prefix+"] "+i_str
        print(o)

    def success(self, i_str):
        o = "SUCCESS ["+self.prefix+"] "+i_str
        print(o)

    def error(self, i_str):
        o = "ERROR ["+self.prefix+"] "+i_str
        print(o)
        raise Exception(o)

    def success_hline(self):
        o = "SUCCESS ["+self.prefix+"] "
        print(o+("*"*(80-len(o))))

    def hline(self):
        print("INFO ["+self.prefix+"] "+"*"*60)



if __name__ == "__main__":

    p = InfoError("test")

