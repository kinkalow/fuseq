import time

class Timer:
    def __init__(self, name):
        self.name = name
        self._start = None
        self._end = None

    def start(self):
        self._start = time.time()

    def end(self):
        self._end = time.time()

    def print(self):
        if not self._start:
            print('[Error] start method is not called')
            return
        end = self._end if self._end else time.time()
        print(f'{self.name} elapased time: {end-self._start:.3f}[s]')
