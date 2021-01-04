
import os
import sys

ADDR_LOC = "./SAMPLING_INFO/"

locs = [ADDR_LOC] # Locations to check/create

for loc in locs:
    if(os.path.isdir(loc)):
        pass
    else:
        os.mkdir(loc)

print("Done!")
