#!/bin/bash

# Usage:
# 1st arg: seq, 2nd arg: chain ID, 3rd arg: output file
# ./seq_to_seqres.sh SKGPGRGDSPYS A output.pdb

# Function to convert one-letter code to three-letter code
convert_aa() {
    case $1 in
        A) echo "ALA";;
        R) echo "ARG";;
        N) echo "ASN";;
        D) echo "ASP";;
        C) echo "CYS";;
        E) echo "GLU";;
        Q) echo "GLN";;
        G) echo "GLY";;
        H) echo "HIS";;
        I) echo "ILE";;
        L) echo "LEU";;
        K) echo "LYS";;
        M) echo "MET";;
        F) echo "PHE";;
        P) echo "PRO";;
        S) echo "SER";;
        T) echo "THR";;
        W) echo "TRP";;
        Y) echo "TYR";;
        V) echo "VAL";;
        *) echo "UNK";;  # Unknown amino acid
    esac
}

# Function to generate SEQRES lines
generate_seqres() {
    local sequence=$1
    local chain=$2
    local length=${#sequence}
    local record=1

    for (( i=0; i<${length}; i+=13 )); do
        # Extract up to 13 residues
        subset=${sequence:$i:13}
        
        # Start the SEQRES line
        printf "SEQRES %3d %s %4d " $record $chain $length
        
        # Convert and print each residue
        for (( j=0; j<${#subset}; j++ )); do
            aa=${subset:$j:1}
            printf " $(convert_aa $aa)"
        done
        
        # New line
        echo ""
        
        # Increment record number
        ((record++))
    done
}

# Check if a sequence is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <sequence> [chain] [output_file]"
    echo "Example: $0 SKGPGRGDSPYS A output.pdb"
    exit 1
fi

# Get the sequence from command line argument
sequence=$1

# Convert to uppercase
sequence=${sequence^^}

# Set the chain identifier (default to A)
chain=${2:-A}

# Generate SEQRES lines
output=$(generate_seqres "$sequence" "$chain")

# If an output file is specified, save to file
if [ -n "$3" ]; then
    echo "$output" > "$3"
    echo "Output saved to $3"
else
    # Otherwise, print to stdout
    echo "$output"
fi