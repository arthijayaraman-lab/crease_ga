class adaptation_params:
    def __init__(self):
        self.gdmmin=0.005
        self.gdmmax=0.85
        self.pcmin=0.5
        self.pcmax=1
        self.pmmin=0.006
        self.pmmax=0.25
        self.kgdm=1.1
        self.pc=0.6
        self.pm=0.001 
    
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