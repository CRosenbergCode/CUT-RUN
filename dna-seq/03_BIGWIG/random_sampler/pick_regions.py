import sys
import random
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-cs", "--contigsizes", help="File with sizes of the genome contigs")
parser.add_argument("-cn", "--contignames", help="File with list of contig names")
parser.add_argument("-r", "--range", help="Range of bases to sample")
parser.add_argument("-n", help="Number of samples to take")
parser.add_argument("-o", "--output", help="Path to output .bed file (include .bed when you enter this)")

args = parser.parse_args()

# Looks like python automatically casts to a long as needed? Gross, but
# it's helpful in this case
pickable_bases = 0

contigs = open(args.contigsizes).readlines()
contig_names = open(args.contignames).readlines()
toremove = []
total_sizes = [0]
for i in range(0, len(contigs)):
    contigs[i] = int(contigs[i].replace('\n', '')) - int(args.range)
    contig_names[i] = contig_names[i].replace('\n', '')
    if contigs[i] <= 0:
        toremove.append(i)
    else:
        pickable_bases += contigs[i]
        total_sizes.append(total_sizes[len(total_sizes) - 1] + contigs[i])
num_removed_indices = 0
for i in toremove:
    contigs.pop(i - num_removed_indices)
    contig_names.pop(i - num_removed_indices)
    num_removed_indices += 1
total_sizes.pop(0)
# print(total_sizes)

def contig_from_location(location):
    if (location > pickable_bases):
        print("error! Selected location %i is greater than total selectable length!" % location)
        exit(1)
    i = 0
    while location > total_sizes[i]:
        i += 1
    return i

def get_picked_regions():
    picked_regions = []
    for i in range(int(args.n)):
        picked_base = random.randint(1, pickable_bases)
        pick = contig_from_location(picked_base)
        picked_name = contig_names[pick]
        picked_start = contigs[pick] - (total_sizes[pick] - picked_base)
        picked_end = picked_start + int(args.range)
        picked_regions.append((picked_name, picked_start, picked_end))
    return picked_regions

def main():
    # Get -n amount of different regions of size -r
    regions = get_picked_regions()

    # Write each region as a .bed line
    towrite = []
    for i in range(len(regions)):
        towrite.append("%s\t%i\t%i\n" % regions[i])
    with open(args.output, 'w') as f:
        f.writelines(towrite)
    
    # print("Sunflower seeds!")

if __name__ == "__main__":
    main()