import argparse


def train_model(training_file, output_file):
    # Read training sequences
    sequences = []
    with open(training_file, "r") as f:
        for line in f:
            sequence = line.strip()
            sequences.append(sequence)

    # Count occurrences of each nucleotide and dinucleotide
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    dinuc_counts = {nt1 + nt2: 0 for nt1 in "ACGT" for nt2 in "ACGT"}
    for seq in sequences:
        for i in range(len(seq)):
            nt = seq[i]
            counts[nt] += 1
            if i < len(seq) - 1:
                dinuc = seq[i : i + 2]
                dinuc_counts[dinuc] += 1

    # Calculate probabilities
    total_counts = sum(counts.values())
    total_dinuc_counts = sum(dinuc_counts.values())
    probs = {}
    for dinuc, count in dinuc_counts.items():
        nt1, nt2 = dinuc[0], dinuc[1]
        prob = count / counts[nt1]
        probs[dinuc] = f"{count}/{counts[nt1]}"

    # Write model to output file
    with open(output_file, "w") as f:
        for dinuc, prob in probs.items():
            f.write(f"{dinuc}\t{prob}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train a zero-order Markov chain model for ESE binding sites"
    )
    parser.add_argument(
        "-f",
        dest="training_file",
        required=True,
        help="Input file with training sequences",
    )
    parser.add_argument(
        "-o", dest="output_file", required=True, help="Output file with trained model"
    )
    args = parser.parse_args()

    train_model(args.training_file, args.output_file)
