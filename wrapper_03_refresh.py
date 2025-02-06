#!/usr/bin/env python3
"""
A CLI wrapper for setting up, downloading, and processing datasets using coderdata.
This script sets up a directory structure, downloads datasets, processes them into
master tables and split files, and outputs various TSV files for downstream analysis.
"""

from copy import deepcopy
import functools as ft
import logging
import argparse
from os import PathLike
from pathlib import Path, PurePath
from typing import Union, List, Optional
import sys

import coderdata as cd
import pandas as pd

# Configure logger
logger = logging.getLogger(__name__)


def main() -> None:
    """Parse command-line arguments and invoke the specified workflow function."""
    main_parser = argparse.ArgumentParser(add_help=True)

    command_parsers = main_parser.add_subparsers(
        dest="command",
        title="commands",
        required=True,
    )

    # Shared arguments for all sub-commands
    p_shared_args = argparse.ArgumentParser(add_help=False)
    p_shared_args.add_argument(
        '-w', '--work_dir',
        dest='WORKDIR',
        type=_check_folder,
        default=Path.cwd(),
        help="Working directory to use (must exist)"
    )
    p_shared_args.add_argument(
        '--overwrite',
        dest='OVERWRITE',
        action='store_true',
        help="Overwrite existing files/folders if they exist"
    )
    p_shared_args.add_argument(
        '-v', '--verbose',
        dest='LOGLEVEL',
        choices=['warn', 'info', 'debug'],
        default='warn',
        help='Defines verbosity level of logging'
    )

    # Sub-command for setting up the workflow
    p_setup_workflow = command_parsers.add_parser(
        "setup",
        parents=[p_shared_args],
        add_help=True,
        help="Set up the required folder structure for the workflow."
    )
    p_setup_workflow.set_defaults(func=setup_workflow)

    # Sub-command for downloading datasets
    p_download_datasets = command_parsers.add_parser(
        "download",
        parents=[p_shared_args],
        add_help=True,
        help="Download all available datasets."
    )
    p_download_datasets.set_defaults(func=download_datasets)

    # Sub-command for processing datasets
    p_process_datasets = command_parsers.add_parser(
        "process",
        parents=[p_shared_args],
        add_help=True,
        help="Process downloaded datasets to generate master tables and splits."
    )
    p_process_datasets.set_defaults(func=process_datasets)
    p_process_datasets.add_argument(
        '-s', '--split_type', dest="SPLIT_TYPE",
        type=str,
        choices=['mixed-set', 'drug-blind', 'cancer-blind'],
        default='mixed-set',
        help="Type of split to perform on the datasets."
    )
    p_process_datasets.add_argument(
        '-n', '--num_splits', dest='NUM_SPLITS',
        type=int,
        default=10,
        help="Number of splits to generate."
    )
    p_process_datasets.add_argument(
        '-r', '--random_seeds', dest='RANDOM_SEEDS',
        type=_random_seed_list,
        default=None,
        help="Comma-separated list of random seeds, one per split."
    )

    # Sub-command to run the full workflow
    p_all = command_parsers.add_parser(
        "all",
        parents=[p_shared_args],
        add_help=True,
        help="Run the full workflow: setup, download, and process datasets."
    )
    p_all.set_defaults(func=full_workflow)

    if len(sys.argv) == 1:
        main_parser.print_help(sys.stderr)
        sys.exit(0)
    try:
        args = main_parser.parse_args()
    except FileNotFoundError as e:
        sys.exit(e)
    except ValueError as e:
        sys.exit(e)

    # Set logging level based on user choice
    if args.LOGLEVEL == 'info':
        loglevel = logging.INFO
    elif args.LOGLEVEL == 'debug':
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARNING

    logging.basicConfig(
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
        level=loglevel
    )
    args.func(args)


def full_workflow(args) -> None:
    """Execute the full workflow: setup, download, and process datasets."""
    setup_workflow(args)
    download_datasets(args)
    process_datasets(args)


def process_datasets(args) -> None:
    """Process the datasets: load data, generate response file, create splits, and build master tables."""
    if args.RANDOM_SEEDS is not None and len(args.RANDOM_SEEDS) != args.NUM_SPLITS:
        sys.exit(
            "<RANDOM_SEEDS> must contain the same number of random seed values as <NUM_SPLITS>."
        )

    local_path = args.WORKDIR.joinpath('data_in_tmp')

    # Retrieve available dataset information
    data_sets_info = cd.list_datasets(raw=True)

    # Load all datasets into a dictionary keyed by dataset name
    logger.info("Importing datasets...")
    data_sets = {}
    for data_set in data_sets_info.keys():
        data_sets[data_set] = cd.load(name=data_set, local_path=local_path)
    logger.info("Importing datasets... done")

    # -------------------------------------------------------------------
    # Concatenate experiment/response data to create response.tsv
    # -------------------------------------------------------------------
    logger.info("Creating 'response.tsv' ...")
    experiments = []
    logger.debug("Creating list of datasets that contain experiment info ...")
    for data_set in data_sets_info.keys():
        if data_sets[data_set].experiments is not None:
            logger.debug(f"Experiment data found for {data_set}")
            experiment = data_sets[data_set].format(
                data_type='experiments',
                shape='wide',
                metrics=[
                    'fit_auc',
                    'fit_ic50',
                    'fit_r2',
                    'fit_ec50se',
                    'fit_einf',
                    'fit_hs',
                    'aac',
                    'auc',
                    'dss',
                ],
            )
            experiments.append(experiment)
        else:
            logger.debug(f"NO experiment data for {data_set}")

    logger.debug("Concatenating experiment data ...")
    response_data = pd.concat(experiments, axis=0, ignore_index=True)
    response_data.rename(
        columns={'improve_drug_id': 'improve_chem_id'},
        inplace=True,
    )
    response_data['improve_sample_id'] = "SAMPLE_ID_" + response_data['improve_sample_id'].astype(int).astype(str)
    outfile_path = args.WORKDIR.joinpath("data_out", "y_data", "response.tsv")
    response_data.to_csv(
        path_or_buf=outfile_path,
        index=False,
        sep='\t',
    )
    logger.info(f"Drug response data written to '{outfile_path}'")
    response_data['index'] = response_data.index  # Temporary index for split reference

    # -------------------------------------------------------------------
    # Create data splits
    # -------------------------------------------------------------------
    split_data_sets(
        args=args,
        data_sets=data_sets,
        data_sets_info=data_sets_info,
        response_data=response_data
    )

    # -------------------------------------------------------------------
    # Get common / reference gene symbols
    # -------------------------------------------------------------------
    data_gene_names = list(data_sets.values())[0].genes
    data_gene_names = (
        data_gene_names[data_gene_names['other_id_source'] == 'ensembl_gene']
        .drop_duplicates(subset=['entrez_id', 'gene_symbol'], keep='first')
    )
    data_gene_names.rename(
        columns={'other_id': 'ensembl_gene_id'},
        inplace=True,
    )
    data_gene_names.drop(columns=['other_id_source'], inplace=True)
    data_gene_names['entrez_id'] = data_gene_names['entrez_id'].astype(int)

    # -------------------------------------------------------------------
    # Create gene expression master table
    # -------------------------------------------------------------------
    merged_transcriptomics = merge_master_tables(
        args=args,
        data_sets=data_sets,
        data_type='transcriptomics'
    )
    merged_transcriptomics = pd.merge(
        merged_transcriptomics,
        data_gene_names[['entrez_id', 'ensembl_gene_id', 'gene_symbol']],
        how='left',
        on='entrez_id',
    )
    merged_transcriptomics.insert(1, 'gene_symbol', merged_transcriptomics.pop('gene_symbol'))
    merged_transcriptomics.insert(0, 'ensembl_gene_id', merged_transcriptomics.pop('ensembl_gene_id'))
    merged_transcriptomics = merged_transcriptomics[merged_transcriptomics['entrez_id'] != 0]
    outfile_path = args.WORKDIR.joinpath("data_out", "x_data", "cancer_gene_expression.tsv")
    (merged_transcriptomics
     .fillna(0)
     .transpose()
     .to_csv(path_or_buf=outfile_path, sep='\t', header=False)
     )

    # -------------------------------------------------------------------
    # Create copy number master table & discretized table
    # -------------------------------------------------------------------
    merged_copy_number = merge_master_tables(args, data_sets=data_sets, data_type='copy_number')
    merged_copy_number.fillna(1, inplace=True)

    # The following bin edges are currently hard-coded; consider moving these to config if needed.
    discretized_copy_number = merged_copy_number.apply(
        pd.cut,
        bins=[0, 0.5210507, 0.7311832, 1.214125, 1.422233, 2],
        labels=[-2, -1, 0, 1, 2],
        include_lowest=True
    )

    merged_copy_number = pd.merge(
        merged_copy_number,
        data_gene_names[['entrez_id', 'ensembl_gene_id', 'gene_symbol']],
        how='left',
        on='entrez_id',
    )
    merged_copy_number.insert(1, 'ensembl_gene_id', merged_copy_number.pop('ensembl_gene_id'))
    merged_copy_number.insert(1, 'gene_symbol', merged_copy_number.pop('gene_symbol'))
    outfile_path = args.WORKDIR.joinpath("data_out", "x_data", "cancer_copy_number.tsv")
    (merged_copy_number
     .transpose()
     .to_csv(path_or_buf=outfile_path, sep='\t', header=False)
     )

    discretized_copy_number = pd.merge(
        discretized_copy_number,
        data_gene_names[['entrez_id', 'ensembl_gene_id', 'gene_symbol']],
        how='left',
        on='entrez_id',
    )
    discretized_copy_number.insert(1, 'ensembl_gene_id', discretized_copy_number.pop('ensembl_gene_id'))
    discretized_copy_number.insert(1, 'gene_symbol', discretized_copy_number.pop('gene_symbol'))
    outfile_path = args.WORKDIR.joinpath("data_out", "x_data", "cancer_discretized_copy_number.tsv")
    (discretized_copy_number
     .transpose()
     .to_csv(path_or_buf=outfile_path, sep='\t', header=False)
     )

    # -------------------------------------------------------------------
    # Create SMILES table
    # -------------------------------------------------------------------
    dfs_to_merge = {}
    for data_set in data_sets:
        if (data_sets[data_set].experiments is not None and data_sets[data_set].drugs is not None):
            dfs_to_merge[data_set] = deepcopy(data_sets[data_set].drugs)

    concat_drugs = pd.concat(dfs_to_merge.values())
    out_df = concat_drugs[['improve_drug_id', 'canSMILES']].drop_duplicates()
    out_df.rename(columns={'improve_drug_id': 'improve_chem_id'}, inplace=True)
    out_df = out_df.dropna(how='any', axis=0)
    outfile_path = args.WORKDIR.joinpath("data_out", "x_data", "drug_SMILES.tsv")
    out_df.to_csv(path_or_buf=outfile_path, sep='\t', index=False)

    # -------------------------------------------------------------------
    # Create mutation count table
    # -------------------------------------------------------------------
    dfs_to_merge = {}
    for data_set in data_sets:
        if data_sets[data_set].experiments is not None and data_sets[data_set].mutations is not None:
            dfs_to_merge[data_set] = data_sets[data_set].mutations
            dfs_to_merge[data_set]['dataset_origin'] = data_set
    merged_mutations = ft.reduce(
        lambda left_df, right_df: pd.merge(
            left_df[['entrez_id', 'improve_sample_id', 'mutation', 'dataset_origin']],
            right_df[['entrez_id', 'improve_sample_id', 'mutation', 'dataset_origin']],
            on=['entrez_id', 'improve_sample_id', 'mutation', 'dataset_origin'],
            how='outer'
        ),
        dfs_to_merge.values()
    )
    unique_mutations = merged_mutations[['entrez_id', 'improve_sample_id', 'mutation']].drop_duplicates()
    unique_mutations['improve_sample_id'] = 'SAMPLE_ID_' + unique_mutations['improve_sample_id'].astype(str)
    mutation_counts = pd.pivot_table(
        unique_mutations,
        values='mutation',
        index='entrez_id',
        columns='improve_sample_id',
        aggfunc='count'
    )
    mutation_counts.fillna(0, inplace=True)
    mutation_counts = pd.merge(
        mutation_counts,
        data_gene_names[['entrez_id', 'ensembl_gene_id', 'gene_symbol']],
        how='left',
        on='entrez_id',
    )
    mutation_counts.insert(1, 'ensembl_gene_id', mutation_counts.pop('ensembl_gene_id'))
    mutation_counts.insert(1, 'gene_symbol', mutation_counts.pop('gene_symbol'))
    mutation_counts = mutation_counts[mutation_counts['gene_symbol'].notna()]
    outfile_path = args.WORKDIR.joinpath("data_out", "x_data", "cancer_mutation_count.tsv")
    mutation_counts.T.to_csv(path_or_buf=outfile_path, sep='\t', header=False)


def split_data_sets(
        args,
        data_sets: dict,
        data_sets_info: dict,
        response_data: pd.DataFrame
) -> None:
    """
    Generate training, testing, and validation splits for each dataset and
    write the indices to separate text files.
    """
    splits_folder = args.WORKDIR.joinpath('data_out', 'splits')
    split_type = args.SPLIT_TYPE
    ratio = (8, 1, 1)
    stratify_by = None
    random_seeds: Optional[List[Optional[int]]] = args.RANDOM_SEEDS if args.RANDOM_SEEDS is not None else [
                                                                                                              None] * args.NUM_SPLITS

    for data_set in data_sets_info.keys():
        if data_sets[data_set].experiments is not None:
            logger.info(f'Creating splits for {data_set} ...')
            drug_response_rows = (
                data_sets[data_set]
                .experiments[['improve_sample_id', 'improve_drug_id', "time", "study"]]
                .drop_duplicates()
            )
            drug_response_rows.rename(
                columns={'improve_drug_id': 'improve_chem_id'},
                inplace=True,
            )
            drug_response_rows['improve_sample_id'] = "SAMPLE_ID_" + drug_response_rows['improve_sample_id'].astype(
                int).astype(str)
            row_nums = pd.merge(
                response_data,
                drug_response_rows,
                how='inner',
                on=['improve_sample_id', 'improve_chem_id', "time", "study"]
            )
            outfile_path = splits_folder.joinpath(f"{data_set}_all.txt")
            row_nums.to_csv(
                path_or_buf=outfile_path,
                columns=['index'],
                index=False,
                header=False
            )

            splits = {}
            for i in range(args.NUM_SPLITS):
                logger.info(f"Split #{i + 1} of {args.NUM_SPLITS} for {data_set} ...")
                splits[i] = data_sets[data_set].train_test_validate(
                    split_type=split_type,
                    ratio=ratio,
                    stratify_by=stratify_by,
                    random_state=random_seeds[i]
                )
                # Process training split
                train_keys = (
                    splits[i].train.experiments[['improve_sample_id', 'improve_drug_id', "time", "study"]]
                    .drop_duplicates()
                )
                train_keys.rename(
                    columns={'improve_drug_id': 'improve_chem_id'},
                    inplace=True,
                )
                train_keys['improve_sample_id'] = "SAMPLE_ID_" + train_keys['improve_sample_id'].astype(int).astype(str)
                row_nums = pd.merge(
                    response_data,
                    train_keys,
                    how='inner',
                    on=['improve_sample_id', 'improve_chem_id', "time", "study"]
                )
                outfile_path = splits_folder.joinpath(f"{data_set}_split_{i}_train.txt")
                row_nums.to_csv(
                    path_or_buf=outfile_path,
                    columns=['index'],
                    index=False,
                    header=False
                )
                logger.debug(f"Training split written to {outfile_path}")

                # Process testing split
                test_keys = (
                    splits[i].test.experiments[['improve_sample_id', 'improve_drug_id', "time", "study"]]
                    .drop_duplicates()
                )
                test_keys.rename(
                    columns={'improve_drug_id': 'improve_chem_id'},
                    inplace=True,
                )
                test_keys['improve_sample_id'] = "SAMPLE_ID_" + test_keys['improve_sample_id'].astype(int).astype(str)
                row_nums = pd.merge(
                    response_data,
                    test_keys,
                    how='inner',
                    on=['improve_sample_id', 'improve_chem_id', "time", "study"],
                )
                outfile_path = splits_folder.joinpath(f"{data_set}_split_{i}_test.txt")
                row_nums.to_csv(
                    path_or_buf=outfile_path,
                    columns=['index'],
                    index=False,
                    header=False
                )
                logger.debug(f"Testing split written to {outfile_path}")

                # Process validation split
                val_keys = (
                    splits[i].validate.experiments[['improve_sample_id', 'improve_drug_id', "time", "study"]]
                    .drop_duplicates()
                )
                val_keys.rename(
                    columns={'improve_drug_id': 'improve_chem_id'},
                    inplace=True,
                )
                val_keys['improve_sample_id'] = "SAMPLE_ID_" + val_keys['improve_sample_id'].astype(int).astype(str)
                row_nums = pd.merge(
                    response_data,
                    val_keys,
                    how='inner',
                    on=['improve_sample_id', 'improve_chem_id', "time", "study"],
                )
                outfile_path = splits_folder.joinpath(f"{data_set}_split_{i}_val.txt")
                row_nums.to_csv(
                    path_or_buf=outfile_path,
                    columns=['index'],
                    index=False,
                    header=False
                )
                logger.debug(f"Validation split written to {outfile_path}")
        logger.info(f"All splits for {data_set} generated")


def merge_master_tables(args, data_sets: dict, data_type: str = 'transcriptomics') -> pd.DataFrame:
    """
    Merge several data tables (e.g., transcriptomics or copy number) into one master table.

    Parameters
    ----------
    args : Namespace
        Command-line arguments.
    data_sets : dict
        Dictionary of dataset objects.
    data_type : str, optional
        Type of data to merge (default is 'transcriptomics').

    Returns
    -------
    pd.DataFrame
        Merged master table.
    """
    dfs_to_merge = []
    for data_set in data_sets:
        if data_sets[data_set].experiments is not None:
            if data_type in ['transcriptomics', 'copy_number'] and getattr(data_sets[data_set], data_type,
                                                                           None) is not None:
                dfs_to_merge.append(
                    data_sets[data_set]
                    .format(data_type=data_type)
                    .transpose()
                    .add_prefix('SAMPLE_ID_', axis=1)
                )
    merged_data: Optional[pd.DataFrame] = None
    for df in dfs_to_merge:
        if merged_data is None:
            merged_data = deepcopy(df)
        else:
            merged_data = merged_data.merge(
                df,
                on='entrez_id',
                suffixes=('', '__rm'),
                how='outer'
            )
            merged_data.columns = merged_data.columns.astype(str)
            merged_data = merged_data.loc[:, ~merged_data.columns.str.contains('__rm')]
    if merged_data is not None and merged_data.index.dtype != int:
        merged_data.index = merged_data.index.astype(int)
    return merged_data


def download_datasets(args) -> None:
    """
    Download all datasets to the temporary data folder.

    Raises
    ------
    FileExistsError
        If the datasets already exist and overwrite is not allowed.
    """
    local_path = args.WORKDIR.joinpath('data_in_tmp')
    exist_ok = args.OVERWRITE
    try:
        cd.download(name='all', local_path=local_path, exist_ok=exist_ok)
    except FileExistsError:
        sys.exit("Data files already exist. Use '--overwrite' to overwrite.")


def setup_workflow(args) -> None:
    """
    Create the folder structure for the workflow:

    <WORKDIR>/
      ├── data_in_tmp     (for downloaded datasets)
      └── data_out
          ├── splits      (for split files)
          ├── x_data      (for master tables)
          └── y_data      (for drug response data)
    """
    parent = args.WORKDIR
    exist_ok = args.OVERWRITE

    data_in = parent.joinpath('data_in_tmp')
    data_out = parent.joinpath('data_out')
    splits = data_out.joinpath('splits')
    x_data = data_out.joinpath('x_data')
    y_data = data_out.joinpath('y_data')

    try:
        data_in.mkdir(exist_ok=exist_ok)
        data_out.mkdir(exist_ok=exist_ok)
        splits.mkdir(exist_ok=exist_ok)
        x_data.mkdir(exist_ok=exist_ok)
        y_data.mkdir(exist_ok=exist_ok)
    except FileExistsError:
        sys.exit(
            "Some folders already exist. To overwrite contents use the '--overwrite' argument."
        )


def _check_folder(path: Union[str, PathLike, Path]) -> Path:
    """
    Validate that the provided path exists and is a directory.

    Parameters
    ----------
    path : Union[str, PathLike, Path]
        The path to check.

    Returns
    -------
    Path
        The absolute path as a Path object.

    Raises
    ------
    TypeError
        If the path is not a string or Path-like object.
    OSError
        If the path does not exist or is not a directory.
    """
    if not isinstance(path, (str, PathLike, Path)):
        raise TypeError(
            f"'path' must be of type str, PathLike or Path. Supplied argument is of type {type(path)}."
        )
    abs_path = path if isinstance(path, Path) else Path(path)
    abs_path = abs_path.absolute()
    if not abs_path.is_dir():
        raise OSError(f"The defined folder path '{path}' does not exist or is not a folder.")
    return abs_path


def _random_seed_list(seed_list_str: str) -> List[int]:
    """
    Convert a comma-separated string of random seeds into a list of integers.

    Parameters
    ----------
    seed_list_str : str
        Comma-separated random seed values.

    Returns
    -------
    List[int]
        A list of integer random seeds.
    """
    seeds = seed_list_str.split(',')
    return [int(item.strip()) for item in seeds]


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
