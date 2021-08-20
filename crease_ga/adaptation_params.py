import numpy as np
class adaptation_params:
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