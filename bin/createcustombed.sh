#!/bin/bash

################################################################################
# Script Name: createcustombed.sh
# Description: Process a panel TSV file based on the specified tool.
# Usage: ./process_panel_tsv.sh <panel_tsv> <tool> 
# Arguments:
#   - panel_tsv: Path to the panel TSV file.
#   - tool:      Tool to use for processing ('oncodrivefml', 'oncodriveclustl', or 'readsxposition').
# Output: An annotated BED file based on the specified tool and header.
# Author: @FedericaBrando
# Date: 2024-05-08
################################################################################

# Function that processes the panel_tsv and merges the regions. 
process_panel() {
    cut -f1,2,6,7 "$1" |
    uniq |
    awk -F'\t' 'NR>1{print $1, $2, $2, $3}' OFS='\t' | 
    bedtools merge -c 4 -o distinct | 
    awk 'BEGIN { FS = "\t"; OFS="\t" } { if ($2 != $3) { $2 += 1 ; $3 -= 1 } else { $2 = $2 ; $3 = $3 }; print }' | #HACK to handle issue arq5x/bedtools2#359
    eval $3
}


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <panel_tsv> <tool>"
    exit 1
fi

panel_tsv="$1"
tool="$2"

# Exit if wrong tool passed
if [ "$tool" != "oncodrivefml" ] && [ "$tool" != "oncodriveclustl" ] && [ "$tool" != "readsxposition" ]; then
    echo "Error: Unsupported tool. Choose from 'oncodrivefml', 'oncodriveclustl', or 'readsxposition'."
    exit 1
fi

# Exits if the panel TSV file does not exist
if [ ! -f "$panel_tsv" ]; then
    echo "Error: Panel TSV file not found: $panel_tsv"
    exit 1
fi

case "$tool" in
    "oncodrivefml")
        HEADER="CHROMOSOME\tSTART\tEND\tELEMENT\tSEGMENT"
        cmd="sed 's/^chr//g' | awk -F'\t' -v HEADER=\"$HEADER\" 'BEGIN {print HEADER} {OFS = \"\t\"} { print \$1, \$2, \$3, \$4, 1}'"
        ;;
    "oncodriveclustl")
        HEADER="CHROMOSOME\tSTART\tEND\tELEMENT_ID\tSYMBOL"
        cmd="sed 's/^chr//g' | awk -F'\t' -v HEADER=\"$HEADER\" 'BEGIN {print HEADER} {OFS = \"\t\"} { print \$1, \$2, \$3, \$4, \$4}'"
        ;;
    "readsxposition")
        HEADER="CHROMOSOME\tSTART\tEND\tGENE\tSYMBOL"
        cmd="awk -F'\t' -v HEADER=\"$HEADER\" 'BEGIN {print HEADER} {OFS = \"\t\"} { print \$1, \$2, \$3, \$4, \$4}'"
        ;;
esac

process_panel "$panel_tsv" "$HEADER" "$cmd"