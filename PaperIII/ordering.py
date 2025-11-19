#!/home/silvia/PhD/bin/python3

import csv
from datetime import datetime
from collections import defaultdict

# Read the file and capture the header
with open("Successful_observations.csv", "r") as f:
    reader = csv.reader(f)
    header = next(reader)  # Save the header
    rows = [row for row in reader]

# Group by frequency (column 2)
groups = defaultdict(list)
for row in rows:
    freq = row[1].strip()
    groups[freq].append(row)

# Sort each group by timestamp (column 3)
for freq in groups:
    groups[freq].sort(key=lambda r: datetime.strptime(r[2].strip(), "%Y-%m-%d %H:%M:%S.%f"))

# Write the output, starting with the header
with open("Successful_observations_sorted.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)  # Write the header first

    # Write 87.68 group first (if it exists)
    if '87.68' in groups:
        for row in groups['87.68']:
            writer.writerow(row)

    # Write other groups
    for freq in sorted(groups):
        if freq != '87.68':
            for row in groups[freq]:
                writer.writerow(row)