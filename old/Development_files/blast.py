import subprocess

# Set paths
blastp_path = r"C:\Program Files\blast\bin\blastp.exe"  # Change to your actual path
query_file = "query.fasta"
database = r"C:\Users\clayt\swissprot"  # Or full path to your database if not in current dir
output_file = "blast_results.txt"

# Construct command
cmd = [
    blastp_path,
    "-query", query_file,
    "-db", database,
    "-out", output_file,
    # sseq qseqid sseqid pident length mismatch gapopen qstart qend ssart send evalue bitscore
    "-outfmt", "6 pident evalue sseqid stitle",     # Tabular format
    #"-evalue","1e-100"
    "-evalue", "10",
    "-task","blastp-short"
]

print("Running BLASTP...")
try:
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    print("BLASTP finished successfully.")
    if result.stdout:
        print("STDOUT:\n", result.stdout)
    if result.stderr:
        print("STDERR:\n", result.stderr)
except subprocess.CalledProcessError as e:
    print("Error running BLASTP!")
    print("Return code:", e.returncode)
    print("Command:", e.cmd)
    print("Output:", e.output)
    print("Error output:\n", e.stderr)

with open("blast_results.txt") as f:
    for i, line in enumerate(f):
        parts = line.strip().split("\t")
        if len(parts) >= 6:
            sseqid = parts[1]
            description = parts[5]  # stitle
            print(f"Top hit [{i+1}]: {sseqid}")
            print("Full description:", description)

            # Optional: Extract just the protein name before OS=
            if "OS=" in description:
                receptor_name = description.split("OS=")[0].strip()
                print("→ Likely receptor:", receptor_name)
            else:
                print("→ Likely receptor:", description)
            break  # Stop after first hit
