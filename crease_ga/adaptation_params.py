import numpy as np
class adaptation_params:
    '''
    Class for all adaptation parameters needed for a GA run.

    Attributes
    ----------
    gdmmin: float. Default=0.005.
        The minimum acceptable value of gdm, a measurement of diversity within
        a generation (high gdm means low diversity, and vice versa). If gdm of
        the current generation falls below `gdmmin`, `pc` will be multiplied by
        `kgdm` and `pm` will be divided by `kgdm` to reduce diversity.
    gdmmax: float. Default=0.85.
        The maximum acceptable value of gdm, a measurement of diversity within
        a generation (high gdm means low diversity, and vice versa). If gdm of
        the current generation exceeds `gdmmax`, `pc` will be divided by
        `kgdm` and `pm` will be multiplied by `kgdm` to increase diversity.
    pcmin: float. Default=0.1.
        Minimum value of `pc`. `pc` cannot be further adjusted below `pcmin`,
        even if `gdm` is still too high.
    pcmax: float. Default=1.
        Maximum value of `pc`. `pc` cannot be further adjusted above `pcmax`,
        even if `gdm` is still too low.
    pmmin: float. Default=0.006.
        Minimum value of `pm`. `pm` cannot be further adjusted below `pmmin`,
        even if `gdm` is still too low.
    pmmin: float. Default=0.25.
        Maximum value of `pm`. `pm` cannot be further adjusted above `pmmax`,
        even if `gdm` is still too high.
    kgdm: float. Default=1.1.
        Should be > 1. The magnitude of adjustment for `pc` and `pm` in case
        `gdm`
        falls ouside of [ `gdmmin`,`gdmmax` ].
    pc: float. Default=0.6.
        possibility of a crossover action happening on an individual in the
        next generation. `pc` is updated after each
        generation according to `gdm`.
    pm: float. Default=0.001.
        possibillity of a mutation action happening on each gene in an
        individual. `pm` is updated after each generation according to `gdm`.
        
    '''
    def __init__(self,
                 gdmmin=0.005,
                 gdmmax=0.85,
                 pcmin=0.1,
                 pcmax=1,
                 pmmin=0.006,
                 pmmax=0.25,
                 kgdm=1.1,
                 pc=0.6,
                 pm=0.001):
        self.gdmmin=gdmmin
        self.gdmmax=gdmmax
        self.pcmin=pcmin
        self.pcmax=pcmax
        self.pmmin=pmmin
        self.pmmax=pmmax
        self.kgdm=kgdm
        self.pc=pc
        self.pm=pm 
    
    def update(self,gdm):
        '''
        Update `pc` and `pm` according to a gdm value.
        '''
        if (gdm > self.gdmmax):
            self.pm *= self.kgdm
            self.pc = np.true_divide(self.pc,self.kgdm)
        elif (gdm < self.gdmmin):
            self.pm = np.true_divide(self.pm,self.kgdm)
            self.pc *= self.kgdm
        if (self.pm > self.pmmax):
            self.pm = self.pmmax
        if (self.pm < self.pmmin):
            self.pm = self.pmmin
        if (self.pc > self.pcmax):
            self.pc = self.pcmax
        if (self.pc < self.pcmin):
            self.pc = self.pcmin
