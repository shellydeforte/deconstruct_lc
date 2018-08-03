import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from itertools import combinations

from deconstruct_lc.analysis_bc.write_bc_score import BcScore
from deconstruct_lc import read_config