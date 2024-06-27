import os
import subprocess
!pip install bio
from Bio import SeqIO
from Bio.Blast import NCBIWWW

import pandas as pd

def parse_blast_text_output(blast_text_file):
    with open(blast_text_file) as file:
        lines = file.readlines()

    data = []
    query_name = None
    bitscore = -1
    evalue = -1
    capture_next = False

    for i, line in enumerate(lines):
        line = line.strip()

        if line.startswith("Query="):
            if query_name is not None:
                # Append the previous query result before processing the new one
                data.append([query_name, bitscore, evalue])
                bitscore = -1
                evalue = -1

            query_name = line.split("Query=")[1].strip()

        if line.startswith(">") and query_name is not None:
            capture_next = True
            continue

        if capture_next and "Score =" in line:
            score_line = line
            bitscore = float(score_line.split("Score = ")[1].split(" bits")[0])
            evalue = float(score_line.split("Expect = ")[1])

            capture_next = False

    if query_name is not None:
        # Append the last query result
        data.append([query_name, bitscore, evalue])

    # Create a DataFrame from the data
    df = pd.DataFrame(data, columns=["Query", "Bitscore", "E-value"])
    return df

# Example usage
blast_text_file = "blast_results"  # Path to your BLAST text output file
df = parse_blast_text_output(blast_text_file)
df.to_csv('blast_output_yay.csv')