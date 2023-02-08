#! /usr/bin/env python3

import socket
import math
import sys
import os
import multiprocessing
import datetime

from mule.InfoError import *
from mule.JobPlatforms import *

__all__ = ['JobParallelizationDimOptions']

class JobParallelizationDimOptions(InfoError):
    """
    Parallelization information for one dimension
    """


    def __init__(self, dim_name = ''):
        """
        Constructor

        Parameters
        ----------
        dim_name: string
            Information on this dimension for error output
        """

        self.init_phase = True

        InfoError.__init__(self, "JobParallelizationDimOptions."+dim_name)

        # Name of dimension for debugging purpose
        self.dim_name = dim_name

        # Overall number of cores
        self.num_cores = None

        # Number of cores per MPI rank.
                # The cores are just allocated, computations are run on threads
        self.num_cores_per_rank = None

        # Number of threads per MPI rank.
                # The threads are used for the computations.
        self.num_threads_per_rank = None

        # Number of ranks
        self.num_ranks = None

        self.init_phase = False


    def __setattr__(self, name, value):

        if name != 'init_phase':
            if not self.init_phase:
                if not name in self.__dict__:
                    raise Exception("Attribute '"+name+"' does not exist!")

        self.__dict__[name] = value


    def print(self):
        self.hline()
        self.info("dim_name: "+str(self.dim_name))
        self.info("num_cores: "+str(self.num_cores))
        self.info("num_cores_per_rank: "+str(self.num_cores_per_rank))
        self.info("num_threads_per_rank: "+str(self.num_threads_per_rank))
        self.info("num_ranks: "+str(self.num_ranks))



    def setup(self, dim_id = None):
        """
        Setup and infer some values automatically.

        In case of some values impossible to infer, generate error
        """

        if self.dim_name == '':
            if dim_id != None:
                self.dim_name = "DIM"+str(dim_id)
                InfoError.setup(self, "JobParallelizationDimOptions."+self.dim_name)

        if self.num_threads_per_rank == None:
            # Use same number of thraeds as there are cores per rank per default
            self.num_threads_per_rank = self.num_cores_per_rank

        if self.num_cores == None:
            if self.num_cores_per_rank != None and self.num_ranks != None:
                self.num_cores = self.num_cores_per_rank*self.num_ranks

        if self.num_cores == None:
            self.error("num_cores not specified")

        if self.num_cores_per_rank == None:
            self.error("num_cores_per_rank not specified")

        if self.num_threads_per_rank == None:
            self.error("num_threads_per_rank not specified")

        if self.num_ranks == None:
            self.error("num_ranks not specified")
        

if __name__ == "__main__":
    p = JobParallelizationDimOptions("testdim")
    p.num_ranks = 12
    p.num_cores_per_rank = 4
    p.setup()
    p.print()

    p.info("FIN")
