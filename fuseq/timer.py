import time

class Timer:
    def __init__(self, name):
        self.name = name

    def start(self):
        self.start = time.time()

    def end(self):
        self.end = time.time()

    def print(self):
        print(f'{self.name} elapased time: {self.end-self.start:.3f}[s]')
