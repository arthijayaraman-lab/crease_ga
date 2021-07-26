class core_params:
    def __init__(self):
        self.popnumber = 5
        self.generations = 10
        self.nloci = 7
        self.minvalu=()
        self.maxvalu=()
        self.numvars = len(self.minvalu)
        
        self.input_file_path = None
        self.q_bounds = None
        self.output_dir = './'