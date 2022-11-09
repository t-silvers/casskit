from abc import ABC, abstractmethod
from dataclasses import MISSING
from pathlib import Path


class Validator(ABC):

    def __set_name__(self, owner, name):
        self.private_name = '_' + name

    def __get__(self, obj, objtype=None):
        return getattr(obj, self.private_name)

    def __set__(self, obj, value):
        self.validate(value)
        setattr(obj, self.private_name, value)

    @abstractmethod
    def validate(self, value):
        pass

class OneOf(Validator):

    def __init__(self, *options):
        self.options = set(options)

    def validate(self, value):
        if value not in self.options:
            raise ValueError(f"Expected {value!r} to be one of {self.options!r}")

class CacheExists(Validator):

    def __init__(self, default: str = MISSING, *options):
        self.default = default

    def __get__(self, obj, obj_type=None):
        return getattr(obj, self.private_name, self.default)

    def validate(self, cache):
        if not Path(cache).exists():
            raise FileNotFoundError(f"Please download biomart data and save it to {cache}")
