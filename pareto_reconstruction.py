import numpy as np
import json
import pareto


data = json.load(open('network_hp_succ.json'))

h_p = cobra.io.load_json_model('iJO1366.json')
