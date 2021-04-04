import ast
import sys
import numpy

file_to_be_written_1 = open(sys.argv[1] + "_list","r")
file_to_be_written_2 = open(sys.argv[2] + "_list","r")

vcflist = open(sys.argv[1] + "_list", "r")
vcflist2 = open(sys.argv[2] + "_list", "r")
annotation_list_1 = open(sys.argv[3], "w")
annotation_list_2 = open(sys.argv[4], "w")


distribution = []
for line in vcflist:
    a_list = ast.literal_eval(line)
    column5 = a_list[3]
    if column5 == ".":
        pass
    else:
        column5 = float(a_list[3])
        distribution.append(column5)

for line in vcflist2:
    a_list = ast.literal_eval(line)
    column5 = a_list[3]
    if column5 == ".":
        pass
    else:
        column5 = float(a_list[3])
        distribution.append(column5)

sdev = numpy.std(distribution)
mean = numpy.mean(distribution)
cutoff = mean + (3*sdev)



for line in file_to_be_written_1:
    a_list = ast.literal_eval(line)
    column1 = str(a_list[0] + "_" + a_list[1])
    column2 = str(a_list[0])
    column3 = str(a_list[1])
    column4 = str(a_list[2])
    column5 = str(a_list[3])
    column5_for_annotation = str(a_list[3])
    if column5 == ".":
        pass
    else:
        column5 = float(a_list[3])
        if column5 < cutoff:
            to_be_written = str(column1 + "\t" + column2 + "\t" + column3 + "\t" + column4 + "\t" + str(column5) + "\n")
        else:
            column5 = str(cutoff)
            to_be_written = str(column1 + "->outlier;real_value=" + column5_for_annotation + "\t" + column2 + "\t" + column3 + "\t" + column4 + "\t" + str(column5) + "\n")
    annotation_list_1.writelines(to_be_written)
    

for line in file_to_be_written_2:
    a_list = ast.literal_eval(line)
    column1 = str(a_list[0] + "_" + a_list[1])
    column2 = str(a_list[0])
    column3 = str(a_list[1])
    column4 = str(a_list[2])
    column5 = str(a_list[3])
    column5_for_annotation = str(a_list[3])
    if column5 == ".":
        pass
    else:
        column5 = float(a_list[3])
        if column5 < cutoff:
            to_be_written = str(column1 + "\t" + column2 + "\t" + column3 + "\t" + column4 + "\t" + str(column5) + "\n")
        else:
            column5 = str(cutoff)
            to_be_written = str(column1 + "->outlier;real_value=" + column5_for_annotation + "\t" + column2 + "\t" + column3 + "\t" + column4 + "\t" + str(column5) + "\n")
    annotation_list_2.writelines(to_be_written)

