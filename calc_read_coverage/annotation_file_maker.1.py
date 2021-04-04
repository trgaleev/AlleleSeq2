import sys

vcf = open(sys.argv[1], "r")
vcf2 = open(sys.argv[2], "r")

list_of_lists = []
for line in vcf:
    stripped_line = line.strip()
    line_list = stripped_line.split()
    list_of_lists.append(line_list)

new_file = open(sys.argv[1] + "_list", "w")

for i in list_of_lists:
    new_file.writelines(str(i)+"\n")

list_of_lists2 = []
for line in vcf2:
    stripped_line = line.strip()
    line_list = stripped_line.split()
    list_of_lists2.append(line_list)

new_file_2 = open(sys.argv[2] + "_list", "w")

for i in list_of_lists2:
    new_file_2.writelines(str(i)+"\n")
