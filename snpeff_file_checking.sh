#!/bin/bash

# Search for snpEff using `which` or `command -v`
SNPEFF_PATH=$(which snpEff)


# Check if snpEff was found
if [[ -z "$SNPEFF_PATH" ]]; then
    echo "Error: snpEff not found in the system's PATH."
    exit 1
else
    echo "snpEff version : $SNPEFF_VERSION"
    echo "snpEff found at: $SNPEFF_PATH"
fi


# Define the path to snpEff data directory
SNPEFF_DIR=$(dirname "$SNPEFF_PATH")
# Extract the conda environment name from the snpEff path


BASE_DIR=$(dirname $(dirname "$SNPEFF_PATH"))

echo "Base directory: $BASE_DIR"


MYCOBACTERIUM_DIR="$BASE_DIR/share/snpeff-5.0-1/data/Mycobacterium_tuberculosis_h37rv"

# Check if the folder exists
if [[ -d "$MYCOBACTERIUM_DIR" ]]; then
    echo "Found Mycobacterium_tuberculosis_h37rv folder at: $MYCOBACTERIUM_DIR"
else
    echo "Error: Mycobacterium_tuberculosis_h37rv folder not found inside snpEff data directory."
    exit 1
fi

SNPEFF_CONFIG="$BASE_DIR/share/snpeff-5.0-1/snpEff.config"
echo $SNPEFF_CONFIG
# Check for the Mycobacterium_tuberculosis_h37rv folder
if [[ ! -d "$MYCOBACTERIUM_DIR" ]]; then
    echo "Mycobacterium_tuberculosis_h37rv folder not found. Downloading..."
    snpEff download -v Mycobacterium_tuberculosis_h37Rv

    # Recheck for the folder
    if [[ ! -d "$MYCOBACTERIUM_DIR" ]]; then
        echo "Error: Failed to download Mycobacterium_tuberculosis_h37Rv data."
        exit 1
    fi
    echo "Downloaded Mycobacterium_tuberculosis_h37rv folder."
else
    echo "Mycobacterium_tuberculosis_h37rv folder found at: $MYCOBACTERIUM_DIR"
fi

# Check for genes.gbk file
if [[ ! -f "$MYCOBACTERIUM_DIR/genes.gbk" ]]; then
    echo "genes.gbk not found in $MYCOBACTERIUM_DIR. Copying provided GenBank file..."
    cp "$REF" "$MYCOBACTERIUM_DIR/genes.gbk"
    echo "Copied $REF as genes.gbk."
fi

# Add genome line to snpEff.config if not already present
GENOME_LINE="AL123456.genome = Mycobacterium_tuberculosis_h37Rv"
if ! grep -Fxq "$GENOME_LINE" "$SNPEFF_CONFIG"; then
    echo "Adding genome line to snpEff.config..."
    sed -i "22i$GENOME_LINE" "$SNPEFF_CONFIG"
    echo "Added genome line to snpEff.config."
else
    echo "Genome line already exists in snpEff.config."
fi

# Build the genome
echo "Building genome..."
snpEff build -c "$SNPEFF_CONFIG" -genbank -v Mycobacterium_tuberculosis_h37rv
if [[ $? -eq 0 ]]; then
    echo "Genome built successfully."
else
    echo "Error: Failed to build the genome."
    exit 1
fi
