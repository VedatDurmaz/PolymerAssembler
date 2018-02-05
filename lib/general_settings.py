#!/usr/bin/env python

import numpy as np


############## --------------------- settings ------------------------------

### config files
usr_cfg_fn = "user_config.yml"
pol_cfg_fn = "lib/polymer_config.yml"

### set label separator for assembled labels
labelSep = '-'

### 3d unity matrix
unity3 = np.array( ((1,0,0),(0,1,0),(0,0,1)) )

### list of values representing true/false
listTrue  = ['true', '1', 't', 'y', 'yes', 'ja', 'wahr', 'yeah']
listFalse = ['false', '0', 'f', 'n', 'no', 'nein', 'falsch']

### list of accepted graphical output formats
listFigure  = ['.eps', '.pdf', '.png']



