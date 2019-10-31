import sys
from collections import defaultdict
import math

input_file = sys.argv[1]
output_file = sys.argv[2]
title = sys.argv[3]
with open(input_file) as file:
    with open(output_file, "a") as out:
        out.write(title)
        total = 0.0

        # build out mapping of taxon id to quantity
        tax_ids = defaultdict(lambda: 0.0)
        for row in file:
            label, sequence = tuple(row.split('\t'))
            tax_ids[label] += 1.0
            total += 1.0

        # write number of unique taxons to file.
        out.write("Number of taxons: " + str(len(tax_ids)))

        # calculate simpson index then write to file.
        numerator = sum([n * (n - 1) for n in tax_ids.values()])
        denominator = total * (total - 1)
        simpson = 1 - (float(numerator) / float(denominator))
        out.write("Simpson Diversity Index: " + str(simpson))

        # calculate shannon index then write to file.
        shannon = -1 * sum((n/total) * math.log(n/total) for n in tax_ids.values())
        out.write("Shannon Diversity Index: " + str(simpson))
