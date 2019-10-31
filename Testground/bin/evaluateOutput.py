import sys
from collections import defaultdict

input_file = sys.argv[1]
output_file = sys.argv[2]
title = sys.argv[3]
with open(input_file) as file:
    with open(output_file, "a") as out:
        out.write(title)
        total = 0

        # build out mapping of taxon id to quantity
        tax_ids = defaultdict(lambda: 0)
        for row in file:
            label, sequence = tuple(row.split('\t'))
            tax_ids[label] += 1
            total += 1

        # write number of unique taxons to file.
        out.write("Number of taxons: " + str(len(tax_ids)))

        # calculate simpson index then write to file.
        numerator = sum([n * (n - 1) for n in tax_ids.values()])
        denominator = total * (total - 1)
        simpson = 1 - (float(numerator) / float(denominator))
        out.write("Simpson Index: " + str(simpson))
