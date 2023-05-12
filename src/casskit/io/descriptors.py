from abc import ABC, abstractmethod
from dataclasses import MISSING
from pathlib import Path


class Validator(ABC):
    
    @abstractmethod
    def validate(self, value):
        pass

class OneOf(Validator):

    def __init__(self, *options):
        self.options = set(options)

    def validate(self, value):
        if value not in self.options:
            raise ValueError(f"Expected {value!r} to be one of {self.options!r}")