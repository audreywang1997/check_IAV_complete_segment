import sys
from Bio import SeqIO
from collections import defaultdict

# Function to check segments for each isolate ID
def check_segments(fasta_file):
    # Define the expected segments
    expected_segments = {"NA", "HA", "NP", "NS", "MP", "PA", "PB1", "PB2"}
    
    # Dictionary to store segments for each isolate ID
    isolate_segments = defaultdict(set)

    # Parse the FASTA file
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            header = record.description
            parts = header.split("|")
            if len(parts) < 6:
                print(f"Warning: Unexpected header format: {header}")
                continue

            isolate_id = parts[1].strip()
            segment = parts[4].strip()

            if segment in expected_segments:
                isolate_segments[isolate_id].add(segment)
            else:
                print(f"Warning: Unexpected segment in header: {header}")

    # Check each isolate ID for missing or extra segments
    results = []
    for isolate_id, segments in isolate_segments.items():
        if len(segments) != 8:
            missing_segments = expected_segments - segments
            results.append((isolate_id, segments, missing_segments))

    return results

# Function to write results to an output file
def write_results(output_file, results):
    with open(output_file, "w") as file:
        for isolate_id, segments, missing_segments in results:
            file.write(f"Isolate ID: {isolate_id}\n")
            file.write(f"Segments present: {', '.join(segments)}\n")
            if missing_segments:
                file.write(f"Missing segments: {', '.join(missing_segments)}\n")
            file.write("\n")

# Main function
def main():
    if len(sys.argv) != 2:
        print("Usage: python check_segments.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = "isolate_segments_report.txt"

    results = check_segments(fasta_file)
    write_results(output_file, results)
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()

