from glob import glob
import shutil
import sys

all_files = glob("RR-SampEn-data/*")
for idx in range(all_files.__len__()):
    new_name = "RR-SampEn-data/subject" + str(idx) + ".mat"
    print(new_name)
    shutil.copy(all_files[idx], new_name)