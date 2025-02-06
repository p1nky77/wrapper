
# IMPROVE Pipeline Data Processing Script

## Overview

This script is a command-line interface (CLI) wrapper that sets up, downloads, and processes datasets for the IMPROVE pipeline using the `coderdata` package. It performs the following main tasks:

- **Setup**: Creates the necessary folder structure.
- **Download**: Downloads all required datasets to a temporary working directory.
- **Process**: Processes the datasets by:
  - Concatenating experiment data into a unified response file.
  - Generating training, testing, and validation splits.
  - Merging transcriptomics and copy number data to form master tables.
  - Generating additional files such as SMILES tables and mutation count tables.

The script utilizes the `tqdm` library to display progress bars in an optimized fashion (using parameters like `miniters` and `tqdm.write()`) and leverages multithreading for I/O-bound dataset loading via `ThreadPoolExecutor`.

---

## Requirements

- **Python**: Version 3.7 or later is recommended.
- **Python Packages**:
  - `coderdata` (ensure you have the correct and trusted version)
  - `pandas`
  - `tqdm`
  - `argparse` (part of the standard library)
- Make sure your working directory has the appropriate permissions for file creation and modification.

---

## Installation

Install the required packages using pip:

```bash
pip install coderdata pandas tqdm
```

---

## Folder Structure

When executed, the script creates a directory structure in the specified working directory (`WORKDIR`):

```
<WORKDIR>/
 ├── data_in_tmp    # Contains downloaded datasets.
 └── data_out
     ├── splits     # Contains split files for each dataset.
     ├── x_data     # Contains combined master tables (e.g., gene expression, copy number).
     └── y_data     # Contains drug response data (response.tsv).
```

---

## Usage

Run the script from the command line with one of the following commands:

```bash
python script_name.py <command> [options]
```

### Available Commands

- **`setup`**  
  Creates the necessary folder structure.

- **`download`**  
  Downloads all required datasets into the `data_in_tmp` folder.

- **`process`**  
  Processes the downloaded datasets. This command:
  - Concatenates experiment data and creates a response file.
  - Generates splits (training, testing, and validation).
  - Merges transcriptomics and copy number data.
  - Creates additional files (SMILES table and mutation count table).

- **`all`**  
  Runs the full workflow (setup, download, and process).

### Command-Line Options

For all commands, the following options are available:

- `-w, --work_dir`:  
  **Description:** Working directory (must exist).  
  **Default:** Current working directory.

- `--overwrite`:  
  **Description:** Overwrite existing files and folders if necessary.

- `-v, --verbose`:  
  **Description:** Sets the logging level.  
  **Choices:** `warn`, `info`, `debug`  
  **Default:** `warn`

Additional options for the `process` command:

- `-s, --split_type`:  
  **Description:** Type of split to generate.  
  **Choices:** `mixed-set`, `drug-blind`, `cancer-blind`  
  **Default:** `mixed-set`

- `-n, --num_splits`:  
  **Description:** Number of splits to generate.  
  **Default:** `10`

- `-r, --random_seeds`:  
  **Description:** Comma-separated list of random seeds. Must have the same length as `NUM_SPLITS`.  
  **Default:** If omitted, random seeds are generated.

### Example

To run the full workflow with informational logging:

```bash
python script_name.py all --verbose info --work_dir /path/to/working_directory
```

---

## Detailed Functionality

### 1. Setup Workflow (`setup` command)
Creates the following directory structure within the specified working directory:
- `data_in_tmp`: For storing downloaded datasets.
- `data_out`: The output directory, containing:
  - `splits`: Generated split files.
  - `x_data`: Master tables for input features (e.g., gene expression, copy number).
  - `y_data`: Drug response data.

### 2. Download Datasets (`download` command)
Uses the `coderdata` package to download all required datasets into the `data_in_tmp` folder. The script checks for existing files and, unless the `--overwrite` flag is set, exits if files already exist.

### 3. Process Datasets (`process` command)
Performs several processing steps:
- **Dataset Import**:  
  Uses a `ThreadPoolExecutor` to load datasets concurrently (optimized for I/O-bound tasks) with progress updates via `tqdm`.
  
- **Experiment Data Processing**:  
  Concatenates experiment data from multiple datasets into a single DataFrame and exports it as `response.tsv` in `data_out/y_data`.
  
- **Split Generation**:  
  Creates training, testing, and validation splits for each dataset. Progress for each dataset and each split is reported using `tqdm` (with `miniters` to reduce update overhead). Messages are output via `tqdm.write()` to avoid interrupting the progress bar.
  
- **Master Table Generation**:  
  Merges and processes transcriptomics and copy number data into master tables. Discretizes copy number data using predefined bin edges and labels.
  
- **Additional Data Files**:  
  Generates a SMILES table from drug data and a mutation count table from mutation data.

### 4. Full Workflow (`all` command)
Executes the entire process (setup, download, process) sequentially.

---

## Logging and Progress Reporting

- **Logging**:  
  The script uses Python’s `logging` module. The logging level is controlled via the `-v/--verbose` option. Log messages are used for error reporting and general status updates.
  
- **Progress Reporting**:  
  The script employs `tqdm` for progress bars. It uses:
  - `miniters` to update the progress less frequently for efficiency.
  - `tqdm.write()` to output messages without disturbing the progress bar.
  - Separate progress bars for dataset loading, experiment data processing, SMILES processing, mutation processing, and splits generation.

---

## Concurrency and Performance

- The dataset import step uses a `ThreadPoolExecutor` to load datasets concurrently, reducing the waiting time for I/O-bound operations.
- `tqdm` is configured to reduce the performance overhead by limiting update frequency (using the `miniters` parameter).

---

## Security Considerations

- **Input Validation**:  
  The script checks that the working directory exists and is a directory.
- **File Path Construction**:  
  Uses `Path.joinpath()` to safely construct file paths, reducing the risk of path injection.
- **Dependencies**:  
  Ensure all third-party packages are installed from trusted sources and are kept up-to-date.
- **Error Handling**:  
  The script logs errors and exits gracefully. Sensitive information should be reviewed before deploying in a production environment.

---

## Customization and Extensibility

- **Adjustable Parameters**:  
  You can modify parameters such as the number of splits, split type, and thread pool size (`max_workers`) to suit your environment.
- **Modular Functions**:  
  The script’s functions (setup, download, process, and split creation) are designed to be modular, so additional processing steps can be integrated as needed.

---

## Conclusion

This script is designed to be robust, efficient, and secure for processing datasets in the IMPROVE pipeline. It incorporates concurrent dataset loading, optimized progress reporting with `tqdm`, and detailed logging. For further customization or troubleshooting, review the inline code documentation and adjust configuration parameters as needed.
