import os
import sys



if __name__ == "__main__": 
    try: 
        file = sys.argv[1]
    except:
        print("Specify file. The product will be stored in the same directory as the file provided.")
        exit()
        
    with open(file, "r") as input_gtf: 
        gtf_lines = input_gtf.readlines()
        
    with open(file.replace(".gtf", ".bed"), "w") as output_bed: 
        for line in gtf_lines: 
            sl = line.split("\t")
            cl_name = sl[-1].split(';')[0].split()[-1].strip('"')
            output_bed.write(f"{sl[0]}\t{int(sl[3])-1}\t{sl[4]}\t{cl_name}\t0\t{sl[6]}\n")
