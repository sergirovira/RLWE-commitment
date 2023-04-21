import os, sys
import csv
import math
from sizes import multiplicative
# import subprocess
# import glob

# parameters are read from filename
filename = "generatedParams.csv"
# if th size of the multiplicativeZKPoK is greater than maxSize that set of parameters is skipped
maxSize = 1024*1024*1024*4

def main():

    args = sys.argv[1:]

    it = 1

    commit = ['keygen','commitment','verifier']
    NIZKP = ['opening','linear','multiplicative']
    options = commit + NIZKP

    with open(filename, newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        headers = next(csv_reader)
        table = list(csv_reader)
        number_params = len(table)

    for opt in (options if args[0] == 'all' else commit if args[0] == 'commit' else NIZKP if args[0] == 'NIZKP' else [args[0]]):
        with open("./benchmarks/skipped.txt", "a") as skipped:
            skipped.write("skipping "+opt+ '\n')
        for j in (range(number_params) if args[1] == '-1' else [int(args[1])]):
            path = "./benchmarks/{0:s}-{1:s}-{2:s}".format(table[j][0],table[j][1],table[j][3])
            if not os.path.exists(path):
                os.makedirs(path)
            for i in range(it):
                print("It: ",i+1)
                print("Q: {0:s} ({1:d})".format(table[j][2],j))
                if multiplicative([float(x) for x in table[j]]) < maxSize:
                    res = os.system("./tests {0:s} {1:d}".format(opt,j))
                    if(res == 0):
                        print("\033[91m" + "[FAIL] something went wrong" +"\033[0m\n")
                        with open("./benchmarks/failed.txt", "a") as failed:
                            failed.write(", ".join(table[j][:4])+" FAIL "+opt)
                        # files = glob.glob("./benchmarks/*")
                        # sorted_by_mtime_ascending = sorted(files, key=lambda t: -os.stat(t).st_mtime)
                        # os.remove(sorted_by_mtime_ascending[0])
                        break
                else:
                    print("\033[91m" + "[SKIP] proof size too large" +"\033[0m\n")
                    with open("./benchmarks/skipped.txt", "a") as skipped:
                        skipped.write(", ".join(table[j][:4])+" -> {0:.2f}GB mult. zkp size\n".format(multiplicative([float(x) for x in table[j]])/(1024*1024*1024*8.)))
                    # we still write into the benchmark folder but logging NaN
                    for name in ([opt] if opt in ['keygen','commitment','verifier']  else ["prover"+opt,"verifier"+opt]):
                        with open("{0:s}/{1:s}_{2:s}.csv".format(path,name,table[j][2]) , "a") as file:
                            file.write("NaN,NaN\n")
                    break


main()
