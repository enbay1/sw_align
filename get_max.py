import sys
max = -float("infinity")
with open(sys.argv[1]) as file_in:
    for line in file_in.readlines():
        test = float(line.strip().split(':')[1])
        if test > max:
            max = test
print(max)