import os

all_files = [x for x in os.listdir() if x.endswith(".txt")]

for file in all_files: 
	with open(file, "r") as inputf: 
		inputf_read = inputf.readlines()
	with open(f"../{file.replace('.txt', '_SMEST1.txt')}", "w") as outputf:
		for line in inputf_read:
			if line.split()[0].endswith("1.1"):
				outputf.write(line)


