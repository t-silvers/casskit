import pandas as pd

from .base import BaseOmic


class CopyNumberVariation(BaseOmic):
    def __init__(
        self,
        data: pd.DataFrame,
    ):
        super().__init__(data)

class MessengerRNA(BaseOmic):
    def __init__(
        self,
        data: pd.DataFrame,
    ):
        super().__init__(data)

class Protein(BaseOmic):
    def __init__(
        self,
        data: pd.DataFrame,
    ):
        super().__init__(data)

class SomaticMutation(BaseOmic):
    def __init__(
        self,
        data: pd.DataFrame,
    ):
        super().__init__(data)
