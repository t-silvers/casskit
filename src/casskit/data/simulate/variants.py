import casskit.data.simulate.base as base


class SimVariants(base.SimulationMixin):

    _data = None

    def __init__(
        self,
        N: int = 100,
        p: int = 500,
        maf_method="gauss",
        **kwargs
    ) -> None:
        super().__init__(N, p)
        self.maf_method = maf_method
        self.kwargs = kwargs
        
        # Methods
        self.methods = {"binary": self.binary}

    @property
    def data(self):
        return self.make_data()

    def make_data(self):
        self.methods[self.cn_method](self.kwargs)

    def binary(self, **kwargs):
        mut_freq = kwargs.get("mut_freq", 0.01)
        return self.rng.choice([0, 1],
                               size=(self.N, self.p),
                               p=[1 - mut_freq, mut_freq])