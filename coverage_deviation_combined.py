import sys
from numpy import median

# Input arguments
int_line_file = sys.argv[1]
wheat_parent_1_file = sys.argv[2]
wheat_parent_2_file = sys.argv[3]
window_size = sys.argv[4]
prefix = sys.argv[5]

# Load coverage data
def load_cov(filepath):
    with open(filepath, 'r') as f:
        return [line.strip().split() for line in f]

int_line = load_cov(int_line_file)
wheat_parent_1 = load_cov(wheat_parent_1_file)
wheat_parent_2 = load_cov(wheat_parent_2_file)

# Chromosome type classifier
def is_alien_chr(chrname):
    return chrname.startswith("C")  # Capital C

def is_wheat_chr(chrname):
    return not chrname.startswith("C")

# Median coverage by type
int_line_alien = [float(x[3]) for x in int_line if is_alien_chr(x[0])]
int_line_wheat = [float(x[3]) for x in int_line if is_wheat_chr(x[0])]

T_median = median(int_line_alien) if int_line_alien else 1
wheat_median = median(int_line_wheat) if int_line_wheat else 1

wp1_wheat_median = median([float(x[3]) for x in wheat_parent_1 if is_wheat_chr(x[0])])
wp2_wheat_median = median([float(x[3]) for x in wheat_parent_2 if is_wheat_chr(x[0])])

# Output list
out_list = []

for i, entry in enumerate(int_line):
    chrname, pos, _, il_cov = entry
    il_cov = float(il_cov)

    if is_alien_chr(chrname):
        # For T or alien chromosomes → IL / T_median
        cov_dev = il_cov / T_median if T_median != 0 else 0
        out_list.append([chrname, pos, str(cov_dev)])

    else:
        # For wheat chrs: IL compared to closest wheat parent
        wp1_cov = float(wheat_parent_1[i][3])
        wp2_cov = float(wheat_parent_2[i][3])
        wp1_norm = wp1_cov / wp1_wheat_median if wp1_wheat_median != 0 else 0
        wp2_norm = wp2_cov / wp2_wheat_median if wp2_wheat_median != 0 else 0
        il_norm = il_cov / wheat_median if wheat_median != 0 else 0

        dev1 = abs(1 - (il_norm / wp1_norm)) if wp1_norm != 0 else float('inf')
        dev2 = abs(1 - (il_norm / wp2_norm)) if wp2_norm != 0 else float('inf')

        if dev1 <= dev2 and wp1_norm != 0:
            cov_dev = il_norm / wp1_norm
        elif wp2_norm != 0:
            cov_dev = il_norm / wp2_norm
        else:
            cov_dev = 0

        out_list.append([chrname, pos, str(cov_dev)])

# Write output
output_file = f"{prefix}_cov_dev_{window_size}.tsv"
with open(output_file, 'w') as out:
    out.write("chr\tposition\tcov_dev\n")
    for row in out_list:
        out.write("\t".join(row) + "\n")
