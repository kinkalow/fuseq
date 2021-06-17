import time

class Timer:
    '''Measure the processing time of a method function'''

    def __init__(self, name):
        self.name = name

    def __call__(self, func):

        def wrap(cls, *args, **kwargs):
            start = time.time()
            ret = func(cls, *args, **kwargs)
            et = time.time() - start
            try:
                # Class dependent
                if cls.params.print_time:
                    self.__print(et)
            except AttributeError:
                self.__print(et)
            return ret
        return wrap

    def __print(self, et):
        print(f'[Info] {self.name} elapsed time: {et:.3f}[s]')
