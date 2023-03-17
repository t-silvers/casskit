from abc import ABC


class BaseOmic(ABC):
    def __init__(self, data):
        self.data = data
        
class BaseMultiOmic(ABC):
    pass