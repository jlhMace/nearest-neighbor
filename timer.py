#### Time class for code timing

import time

class Timer:
    # Class variable for all Timer instances
    timers = {}

    def __init__(
            self,
            name = None,
            logger = print
            ):
        self.start_time = None
        self.name = name
        self.logger = logger
        self.total_time = 0

        # Add named timer to dict of timers
        if name:
            self.timers.setdefault(name,0)
        
    def start(self):
        '''Start timer'''
        self.start_time = time.perf_counter()

    def stop(self):
        '''Stop timer'''
        elapsed_time = time.perf_counter() - self.start_time

        #if self.logger:
        #    self.logger(f"    Elapsed time: {elapsed_time}")
        if self.name:
            self.timers[self.name] += elapsed_time
        self.start_time = None
        self.total_time += elapsed_time

    def return_time(self):
        '''Return elapsed time'''
        return self.total_time
