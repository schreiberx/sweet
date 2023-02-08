#! /usr/bin/env python3

import sys
import os
import concurrent.futures
import queue
import multiprocessing
import ctypes


#
# Helper routine for parallelization across ranges
#
# This is based on worker threads which will start the tasks across a specified range
#


class parhelper:
    def __init__(self):
        pass


    def _worker(self, id):
        if self.verbose > 0:
            print("Hello from worker "+str(id))

        work_items = 0
        while True:
            # Get next work id
            with self.counter.get_lock():
                c = self.counter.value
                if c >= self.range_size:
                    c = None
                else:
                    self.counter.value += 1

            if c == None:
                break

            # Call actual function
            self.fun(c)

            if self.use_tqdm:
                # Once finished, update progress bar
                #pbar = self.q.get()
                self.pbar.update()
                #self.q.put(pbar)

            work_items += 1

        if self.verbose > 0:
            print("Worker "+str(id)+" processed "+str(work_items)+" work items")



    def __init__(self, range_size, fun, max_workers=None, use_tqdm=True, verbose=0):
        self.range_size = range_size
        self.fun = fun
        self.use_tqdm = use_tqdm

        self.verbose = verbose
        self.max_workers = max_workers

        if self.max_workers == None:
            self.max_workers = os.cpu_count()

        self.threadpool = concurrent.futures.ThreadPoolExecutor(max_workers)

        # Atomic counter
        self.counter = multiprocessing.Value(ctypes.c_int64, 0)

        # Prepare queue to hand over pbar handler
        #self.q = queue.Queue()

        if use_tqdm:
            import tqdm
            self.pbar = tqdm.tqdm(total=self.range_size)
            #self.q.put(pbar)

        # Enqueue the worker threads which will do the work stealing for us
        with self.threadpool:
            self.threadpool.map(self._worker, range(self.max_workers))

        if use_tqdm:
            #pbar = self.q.get()
            self.pbar.close()



if __name__ == "__main__":

    print("Constructor")

    def fun(i):
        #print(" + work item "+str(i))
        import time
        time.sleep(0.00001)


    n = 100000
    print("Run "+str(n)+" tasks")
    parhelper(n, fun)
    parhelper(n, fun, 4, use_tqdm=True, verbose=1)

