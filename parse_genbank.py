import argparse
from Bio import SeqIO

def parse_genbank(file_path, output_file):
    # Open output file
    with open(output_file, 'w') as out:
        # Parse the GenBank file
        for record in SeqIO.parse(file_path, "genbank"):
            # Loop through each feature in the record
            for feature in record.features:
                feature_type = feature.type
                start = int(feature.location.start)
                end = int(feature.location.end)

                # Try to retrieve the name from 'label', 'note', 'product', then fallback to 'gene'
                feature_name = feature.qualifiers.get('label') or \
                               feature.qualifiers.get('note') or \
                               feature.qualifiers.get('product') or \
                               feature.qualifiers.get('gene') or ['Unnamed feature']

                # Write to the output file in CSV format
                out.write(f"{feature_name[0]},{feature_type},{start},{end}\n")

if __name__ == "__main__":
    # Parse command line arguments for input and output file
    parser = argparse.ArgumentParser(description="Parse GenBank file and output CSV")
    parser.add_argument("input_file", help="Path to the input GenBank file (.gb)")
    parser.add_argument("output_file", help="Path to the output CSV file")
    
    args = parser.parse_args()

    # Call the parsing function with the provided files
    parse_genbank(args.input_file, args.output_file)
