# /genetics/PAREG/perrotn/scripts/combinations_MVMR.py

##########################################################
###### GENERATE ALL EXPOSURES COMBINATIONS FOR MVMR ######
##########################################################

# 1) Initialize
import argparse

# Instantiate the parser
parser = argparse.ArgumentParser(description='Generate all combinations for MVMR')
# 2) Add Arguments

# Required positional argument
parser.add_argument('--file', type=str, help='Filename for MVMR')


# 3) Parse
args = parser.parse_args()

# 4) Access
print("Filename:")
print(args.file)

import numpy as np
import pandas as pd
from itertools import combinations
import re
import os

# TEST
# filename="/storage/genetics_work2/perrotn/PURE/MVMR_consortia/TEST_FILE_MVMR.txt"

filename=args.file

df=pd.read_csv(filename, sep = "\t")

# for i in range(1, len(df.exposures)+1):
#     try:
#         os.remove(re.sub(".txt", "_"+str(i)+".txt", filename))
#     except OSError:
#         pass

[pd.DataFrame(list(combinations(df.exposures, i))).to_csv(path_or_buf=re.sub(".txt", "_"+str(i)+".txt", filename), sep="\t", header=False, index=False) for i in range(1, len(df.exposures)+1)]


