#Author: Lucas Zenichi Terada
#institution: University of Campinas
import pandas as pd
import json

bus = None
branches = None
switches = None
param = None

def loadfiles(self):
    self.bus = pd.read_csv('/home/zenichi/Dropbox/Research/2020 - Fapesp Scholarship/dnotools/example/powerflow/3Bus/bus.csv')
    self.branches = pd.read_csv('/home/zenichi/Dropbox/Research/2020 - Fapesp Scholarship/dnotools/example/powerflow/3Bus/networks.csv')
    self.switches = pd.read_csv('/home/zenichi/Dropbox/Research/2020 - Fapesp Scholarship/dnotools/example/powerflow/3Bus/switches.csv')
    path = '/home/zenichi/Dropbox/Research/2020 - Fapesp Scholarship/dnotools/example/powerflow/3Bus/param.json'
    with open(path) as json_file:
        self.param= json.load(json_file)