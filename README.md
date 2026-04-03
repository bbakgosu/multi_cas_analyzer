# Multi-Cas Analyzer

A Python tool for analyzing CRISPR-Cas editing outcomes from FASTQ files. Automatically detects insertion, deletion, substitution, and complex edit patterns, and generates Excel reports.

## Features

- Automated CRISPR editing pattern analysis from FASTQ files
- Classification of Wild Type (WT), Insertion, Deletion, Substitution, and Complex edits
- Excel report generation with per-edit-type sheets
- Sequence alignment visualization
- Multi-target analysis support

## Requirements

```bash
pip install -r requirements.txt
```

Required packages:
- biopython
- pandas
- openpyxl
- regex

## Usage

### Basic

```bash
python multi_cas_analyzer.py <fastq_folder> <parameter_file> <edit_range>
```

### Arguments

- `fastq_folder`: Path to the folder containing FASTQ files
- `parameter_file`: Parameter file with target information
- `edit_range`: Analysis window around the cleavage site (bp)

### Example

```bash
python multi_cas_analyzer.py ./fastq NovaIscB_VEGFA_wRNA1_param.txt 10
```

## Parameter File Format

Each target is defined by three consecutive lines:

```
Target_Name
REFERENCE_SEQUENCE
SPACER_SEQUENCE
```

Example:
```
NovaIscB_VEGFA_wRNA1
ATCGATCGATCG...
ATCGATCG
```

## Output

The program writes the following files to the `output` folder:

1. **Text report** (`*.fastq.txt`)
   - Per-FASTQ summary statistics
   - WT, Insertion, Deletion, Substitution, and Complex counts and percentages

2. **Excel report** (`*_patterns.xlsx`)
   - Summary sheet: overview of all edit types
   - Individual sheets per edit type with alignment and pattern counts
   - Color-coded formatting per edit category

3. **Summary file** (`summary.txt`)
   - Aggregated results across all FASTQ files

## Analysis Workflow

1. Read FASTQ files
2. Crop sequences using indicator sequence
3. Align against the reference sequence
4. Identify variant patterns within the edit range
5. Classify and count edit types
6. Generate reports

## Edit Type Definitions

| Type | Description |
|------|-------------|
| **WT** | Wild type — no variants detected |
| **Insertion** | Base insertion |
| **Deletion** | Base deletion |
| **Substitution** | Base substitution |
| **Complex** | Combination of multiple edit types |

## Configuration

Adjust the `MINIMUM_FREQUENCY` variable to set the minimum read count threshold (default: 2).

## License

Made by Chanju
