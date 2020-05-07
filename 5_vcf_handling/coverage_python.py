import subprocess
import pandas as pd


direct="{}/{}/vcf_handling/".format(path, sp)

# Function call
def average(direct, name):
    f = open("{}{}".format(direct, name))
    average = 0
    Sum = 0
    row_count = 0
    for line in f:
        dp = line.split()[0]
        n=float(dp)
        Sum += n
        row_count += 1
    average = Sum / row_count
    average = round(average,2)
    summary = "{} {} {}".format(name, row_count, average)
    f.close()
    print(summary)
