#!/usr/bin/env python3.6

import subprocess

weight_and_dont_care_position_combinations = """
|2|2|
|2|3|
|2|4|
|2|5|
|2|6|
|2|7|
|2|8|
|2|9|
|2|10|
|2|11|
|2|12|
|2|13|
|2|14|
|2|15|
|2|16|
|3|3|
|3|4|
|3|5|
|3|6|
|3|7|
|3|8|
|3|9|
|3|10|
|3|11|
|3|12|
|3|13|
|3|14|
|3|15|
|3|16|
|4|3|
|4|4|
|4|5|
|4|6|
|4|7|
|4|8|
|4|9|
|4|10|
|4|11|
|4|12|
|4|13|
|4|14|
|4|15|
|4|16|
|5|3|
|5|4|
|5|5|
|5|6|
|5|7|
|5|8|
|5|9|
|5|10|
|5|11|
|5|12|
|5|13|
|5|14|
|5|15|
|5|16|
|7|3|
|7|4|
|7|5|
|7|6|
|7|7|
|7|8|
|7|9|
|7|10|
|7|11|
|7|12|
|7|13|
|7|14|
|7|15|
|7|16|
"""

rasbhari_path = "/home/robin/code/rasbhari-v1.4.0/rasbhari"
patterns_per_set = ['-m', '5']
background_match_prob = ['-q', '0.05']
match_prob = ['-p', '0.3']
sequence_length = ['-S', '400']
homologue_region_length = ['-H', '20']


for line in weight_and_dont_care_position_combinations.strip().split('\n'):
    if line.startswith("#"):
        continue
    weight, dont_care = [t for t in line.split('|') if t]
    print(f'Generating patterns for weight {weight} and don\'t care {dont_care}')

    subprocess.run(
        [rasbhari_path, '--nosens', *patterns_per_set, *background_match_prob,
         *match_prob, *sequence_length, *homologue_region_length,
         '-w', weight, '-d', dont_care, '--outfile', f'./data/w-{weight}_d-{dont_care}.out'])