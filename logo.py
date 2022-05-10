"""
This script is a wrapper around Roni's motif creation workflow.

It can be used in two ways:
1) Create a single motif and save it as png file.
2) Create a matplotlib/logomaker plot with number of motifs (e.g. all_a2g, a2g_in_repeats) and return it, or save to
the disk too.
"""

import argparse
import inspect
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from logomaker import Logo

def create_fasta_file(
    in_bed_file,
    genome_file,
    surrounding_bases,
    temp_dir,
    fasta_file,
    motif_col_position,
    skip_header,
    stranded_apobec_motif=False
):
    """
    create a fasta file, each line containing a seq of an editing site + $surrounding_bases from each side.

    Save the fasta file to temp_dir.

    TODO: re-implement without relying on external installations of `awk` and `bedtools`

    @param in_bed_file: bed file to create the motif from
    @type in_bed_file: Path
    @param genome_file: genome fasta file
    @type genome_file: Path
    @param surrounding_bases: number of bases from each side of edited positions (e.g. 1 for ADAR, 3 for APOBEC)
    @type surrounding_bases: int
    @param temp_dir: temp_bed_file and fasta_file are saved in it
    @type temp_dir: Path
    @param fasta_file: out path of fasta file created from in_bed_file
    @type fasta_file: Path
    @param motif_col_position: used for editing files of old editing_detection, should be removed
    @type motif_col_position: int
    @param skip_header: a flag indicating the first line of bed_file if a comment (e.g. #chrom start end)
    @type skip_header: bool
    """
    temp_bed_file = Path(temp_dir, f"temp_bed_file.bed")
    with temp_bed_file.open("w") as out_file:
        with in_bed_file.open("r") as in_file:
            for i, line in enumerate(in_file, start=1):
                if i == 1 and skip_header:
                    continue
                line_parts = line.split("\t")
                line_parts[3] = line_parts[3].split(";")[
                    motif_col_position
                ]  # extract wanted motif (ADAR or APOBEC)
                if stranded_apobec_motif:
                    line_parts[5] = "+\n" if line_parts[5] == "-\n" else "-\n"
                new_line = "\t".join(line_parts[x] for x in range(len(line_parts)))
                out_file.write(new_line)
    # the awk command sets the new coordinates
    # if stranded_apobec_motif:
    # awk -v OFS='\t' '{if ($6=="-") {print $1, $2-3, $3+3, $4, $5, "+"} else {print $1, $2-3, $3+3, $4, $5, "-"}}'
    awk_cmd = (
            f"awk -v OFS='\t' '{{print $1, $2-{surrounding_bases}, $3+{surrounding_bases}, $4, $5, $6}}' "
            f"{temp_bed_file}"
        )
    # the bedtools cmd creates a fasta file according to the extended coordinates from the prev awk cmd
    bedtools_cmd = f"bedtools getfasta -s -fi {genome_file} -bed stdin"
    final_cmd = f"{awk_cmd} | {bedtools_cmd} > {fasta_file}"
    subprocess.run(final_cmd, shell=True)


# head -n 1 Data/HyperEditing/ClusterScreening.SE_0.05_0.6_30_0.6_0.1_0.8_0.2/Original/ESuniqS.stranded.A2G.bed |\
# awk -v OFS='\t' '{print $1, $2-3, $3+3, $4, $5, $6}'

# NC_011814.1	2475	2476	TCT	1	-
# NC_029649.1	11866285	11866286	CGG	1	+

# >NC_011814.1:2472-2479(-)
# CTAGAAG
# >NC_029649.1:11866282-11866289(+)
# GACGGCG

# cat Data/HyperEditing/ClusterScreening.SE_0.05_0.6_30_0.6_0.1_0.8_0.2/Original/ESuniqS.stranded.C2T.bed |\
# awk -v OFS='\t' '{print $1, $2-3, $3+3, $4, $5, if ($6)}' |\
# bedtools getfasta -s -fi Data/Annotations/genome.fasta -bed stdin |\
# less -N

# cat Data/HyperEditing/ClusterScreening.SE_0.05_0.6_30_0.6_0.1_0.8_0.2/Original/ESuniqS.stranded.C2T.bed |\
# awk -v OFS='\t' '{print $1, $2-3, $3+3, $4, $5, $6}' |\
# bedtools getfasta -s -fi Data/Annotations/genome.fasta -bed stdin |\
# less -N


def create_frequency_file(temp_dir, fasta_file):
    final_cmd = (
        f"python2.7 /home/alu/fulther/Scripts/scripts_python/General/LogoFreqAnalyzer.py -o {temp_dir} -i "
        f"{fasta_file}"
    )
    subprocess.run(final_cmd, shell=True)


def plot_single_logo(temp_dir):
    """Plot a single logo with SequenceMotifFreq.tab in temp_dir."""
    final_cmd = (
        f"/usr/bin/Rscript /home/alu/fulther/Scripts/scripts_R/General/ggseqLOGO.R -i "
        f"{Path(temp_dir, 'SequenceMotifFreq.tab')} -o {temp_dir} --png"
    )
    subprocess.run(final_cmd, shell=True)


def clean_up(temp_dir, output_file, keep_temp_files):
    logo_file = Path(temp_dir, "BaseFreqLOGO.png")
    cp_cmd = f"cp -T {logo_file} {output_file}"
    subprocess.run(cp_cmd, shell=True)
    if not keep_temp_files:
        delete_folder_with_files(temp_dir)


# def plot_multiple_logos(temp_dirs, titles, main_title=None, surrounding_bases=1):
#     """
#     Plot a figure with multiple motifs of ADAR out of the frequency files in temp_dirs.
#
#     # todo: use modulo s.t. 4 motifs will produce a 2x2 fig, 5 motifs a 3x2 fig where the bottom-right is missing, etc.
#
#     @param temp_dirs: dirs with SequenceMotifFreq.tab in them, each following the template
#     /some_path/bed_file_name.temp_logo_dir
#     @type temp_dirs: list
#     @param titles: titles of sub plots
#     @type titles: list
#     @param main_title: main title of the plot, default=""
#     @type main_title: str or None
#     @return: fig
#     @rtype: matplotlib.figure.Figure
#     """
#     # define titles
#     titles = (
#         titles if titles else [temp_dir.name.split(".")[0] for temp_dir in temp_dirs]
#     )
#     # if not titles:
#     #     titles = [temp_dir.name.split(".")[0] for temp_dir in temp_dirs]
#     # define dfs with frequencies
#     frequency_dfs_paths = [
#         Path(temp_dir, "SequenceMotifFreq.tab") for temp_dir in temp_dirs
#     ]
#     frequency_dfs = [
#         pd.read_csv(frequency_df_path, sep="\t").iloc[:, 1:]
#         for frequency_df_path in frequency_dfs_paths
#     ]
#     # define fig, sub plots and titles
#     fig_width = 4 * len(temp_dirs)
#     fig_size = (fig_width, 3)
#     num_of_plots = len(temp_dirs)
#     fig, axes = plt.subplots(
#         nrows=1,
#         ncols=num_of_plots,
#         figsize=fig_size,
#         gridspec_kw={"wspace": 0.5},
#         squeeze=False,
#     )
#     if main_title:
#         fig.subplots_adjust(top=0.8)
#         # fig.suptitle(main_title, y=1.0, fontsize="x-large")
#         fig.suptitle(main_title, fontsize="x-large")
#     plots = []
#     for i, (title, frequency_df) in enumerate(zip(titles, frequency_dfs)):
#         # ax = axes[i]
#         ax = axes[0][i]
#         ax.set_title(title, pad=15)
#         ax.set_title(title)
#         plot = Logo(df=frequency_df, ax=ax, vpad=0.02)
#         plots.append(plot)
#     # style plots
#     for plot in plots:
#         # style plots using Logo methods
#         plot.style_spines(visible=False)
#         plot.style_spines(spines=("left", "bottom"), visible=True)
#         # style plots using matplotlib's Axes methods
#         plot.ax.set_ylabel("Probability", labelpad=5)
#         plot.ax.set_xlabel("Position", labelpad=5)
#
#         # plot.ax.set_xticks([0, 1, 2])
#         # plot.ax.set_xticklabels(["-1", "0", "1"])
#         plot.ax.set_xticks(list(range(surrounding_bases * 2 + 1)))
#         plot.ax.set_xticklabels(
#             list(range(-surrounding_bases, 0, 1))
#             + [0]
#             + list(range(1, surrounding_bases + 1))
#         )
#
#     # todo check if this helps to show x-axis title
#     _, top_ylim = plt.ylim()
#     plt.ylim(0, top_ylim)
#
#     return fig


# def multiple_logos_heatmap(
#     bed_files: list[Path],
#     temp_dirs: list[Path],
#     genome_file: Path,
#     surrounding_bases=1,
#     keep_temp=True,
#     motif_col_position: int = 0,
#     skip_header: bool = False,
#     titles: list[str] = None,
#     stranded_apobec_motif:bool = False
# ):
#     # (1) set up paths
#     for temp_dir in temp_dirs:
#         temp_dir.mkdir(exist_ok=True)
#     genome_file = genome_file.absolute()
#     # (2) create a fasta file for each bed_file in bed_files, where each line contains a seq of
#     # surrounding_bases-edited_base-surrounding_bases for each edited_base in the bed_file, and save it to temp_dir
#     fasta_files = [
#         Path(temp_dir, "positions_from_sides.fasta") for temp_dir in temp_dirs
#     ]
#     for bed_file, fasta_file, temp_dir in zip(bed_files, fasta_files, temp_dirs):
#         create_fasta_file(
#             bed_file,
#             genome_file,
#             surrounding_bases,
#             temp_dir,
#             fasta_file,
#             motif_col_position,
#             skip_header,
#             stranded_apobec_motif
#         )
#     # (3) create a frequency file out of each fasta_file in fasta_files, and save it to that fasta_file's corresponding
#     # temp_dir
#     for fasta_file, temp_dir in zip(fasta_files, temp_dirs):
#         create_frequency_file(temp_dir, fasta_file)
#     frequency_files = [
#         Path(temp_dir, "SequenceMotifFreq.tab") for temp_dir in temp_dirs
#     ]
#
#     # define titles
#     titles = (
#         titles
#         if titles
#         else [temp_dir.name.split("-")[0] for temp_dir in temp_dirs]
#     )
#     # return frequency_files, titles
#
#     frequency_dfs = []
#     for title, frequency_file in zip(titles, frequency_files):
#         df = pd.read_csv(frequency_file, sep="\t")
#         df.insert(1, "Sample", title)
#         frequency_dfs.append(df)
#     return frequency_dfs


# def multiple_logos(
#     bed_files,
#     genome_file,
#     titles=None,
#     output_file="",
#     surrounding_bases=1,
#     motif_col_position=0,
#     keep_temp_files=False,
#     skip_header=False,
#     main_title=None,
#     stranded_apobec_motif=False
# ):
#     """
#     Create multiple logos out of bed_files and plot them together, side by side.
#
#     This is one of the two main functions of the module.
#     It creates a matplotlib/logomaker plot with number of motifs (e.g. all_a2g, a2g_in_repeats) and returns it.
#     Optionally, by providing a nonempty string for output_file, the plot can be saved to disk too.
#
#     @param bed_files: bed files
#     @type bed_files: list[Path]
#     @param genome_file: genome fasta file
#     @type genome_file: Path
#     @param titles: titles of each motif in the plot - inferred from one-before-last suffix (e.g 'ESuniqS.A2G.bed'
#     will result in an 'A2G' subtitle) if subtitles is None
#     @type titles: list[str]
#     @param output_file: save the plot to output_file if output_file is provided
#     @type output_file: Path
#     @param surrounding_bases: number of bases from each side of edited positions (e.g. 1 for ADAR, 3 for APOBEC)
#     @type surrounding_bases: int
#     @param motif_col_position: used for editing files of old editing_detection, should be removed
#     @type motif_col_position: int
#     @param keep_temp_files: a flag indicating whether to keep temp files, only works in conjunction with output_file
#     @type keep_temp_files: bool
#     @param skip_header: a flag indicating the first line of bed_file if a comment (e.g. #chrom start end)
#     @type skip_header: bool
#     @param main_title: main title of the plot, default=""
#     @type main_title: str or None
#     @return: fig
#     @rtype: matplotlib.figure.Figure
#     """
#     # (1) set paths
#     temp_dirs_base_path = (
#         output_file.parent if isinstance(output_file, Path) else bed_files[0].parent
#     )
#     temp_dirs = [
#         Path(temp_dirs_base_path, f"{bed_file.stem.split('.')[-1]}.temp_logo_dir")
#         for bed_file in bed_files
#     ]
#     for temp_dir in temp_dirs:
#         temp_dir.mkdir(exist_ok=True)
#     # (2) create a fasta file for each bed_file in bed_files, where each line contains a seq of
#     # surrounding_bases-edited_base-surrounding_bases for each edited_base in the bed_file, and save it to temp_dir
#     # fasta_files = [Path(temp_dir, "positions_from_sides.bed") for temp_dir in temp_dirs] TODO check
#     fasta_files = [
#         Path(temp_dir, "positions_from_sides.fasta") for temp_dir in temp_dirs
#     ]
#     for bed_file, fasta_file, temp_dir in zip(bed_files, fasta_files, temp_dirs):
#         create_fasta_file(
#             bed_file,
#             genome_file,
#             surrounding_bases,
#             temp_dir,
#             fasta_file,
#             motif_col_position,
#             skip_header,
#             stranded_apobec_motif
#         )
#     # (3) create a frequency file out of each fasta_file in fasta_files, and save it to that fasta_file's corresponding
#     # temp_dir
#     for fasta_file, temp_dir in zip(fasta_files, temp_dirs):
#         create_frequency_file(temp_dir, fasta_file)
#     # (4) plot the motif according to the frequency file, and save it to temp_dir
#     fig = plot_multiple_logos(
#         temp_dirs, titles, main_title=main_title, surrounding_bases=surrounding_bases
#     )
#     # (5) delete temp_dir if needed, save the fig (if needed) and return it in any case
#     if not keep_temp_files or not output_file:
#         delete_folders_with_files(temp_dirs)
#     scale_fig_width_inches_by = 1.25
#     scale_fig_height_inches_by = 1.25
#     fig = plt.gcf()
#     width, height = fig.get_size_inches()
#     width, height = (
#         width * scale_fig_width_inches_by,
#         height * scale_fig_height_inches_by,
#     )
#     fig.set_size_inches(width, height)
#     fig.tight_layout()
#     if output_file:
#         fig.savefig(output_file)
#     return fig


def single_logo(
    bed_file,
    genome_file,
    output_file,
    surrounding_bases,
    motif_col_position,
    keep_temp_files,
    skip_header,
):
    """
    Create a single logo out of bed_file.

    This is one of the two main functions of the module.
    It is mainly a wrapper around Roni's scripts, and can be run either by importing it or through the command
    line.
    The logo is saved as a .png file.

    @param bed_file: tested bed file
    @type bed_file: Path
    @param genome_file: genome fasta file
    @type genome_file: Path
    @param output_file: output path of plotted logo
    @type output_file: Path
    @param surrounding_bases: number of bases from each side of edited positions (e.g. 1 for ADAR, 3 for APOBEC)
    @type surrounding_bases: int
    @param motif_col_position: used for editing files of old editing_detection, should be removed
    @type motif_col_position: int
    @param keep_temp_files: a flag indicating whether to keep temp files
    @type keep_temp_files: bool
    @param skip_header: a flag indicating the first line of bed_file if a comment (e.g. #chrom start end)
    @type skip_header: bool
    """
    # (1) set paths
    temp_dir = Path(output_file.parent, f"{os.getpid()}_temp_logo_dir")
    temp_dir.mkdir()
    # (2) create a fasta file, where each line contains a seq of a site from the bed file surrounded by
    # $surrounding_bases from each side, and save it to temp_dir
    fasta_file = Path(temp_dir, "positions_from_sides.bed")
    create_fasta_file(
        bed_file,
        genome_file,
        surrounding_bases,
        temp_dir,
        fasta_file,
        motif_col_position,
        skip_header,
    )
    # (3) create a frequency file out of the fasta file, and save it to temp_dir
    create_frequency_file(temp_dir, fasta_file)
    # (4) plot the motif according to the frequency file, and save it to temp_dir
    plot_single_logo(temp_dir)
    # (5) copy the plot from temp_dir to the predefined final output_file, and delete temp_dir if not keep_temp_files
    clean_up(temp_dir, output_file, keep_temp_files)

def delete_folders_with_files(folders):
    """
    Delete each folder and its files.

    @param folders: folders with files
    @type folders: iterable
    """
    for folder in folders:
        delete_folder_with_files(folder)


def delete_folder_with_files(folder):
    """
    Delete folder and its files.

    @param folder: folder with files to delete
    @type folder: str or Path
    """
    subprocess.run(f"rm -rf {folder}", shell=True)


if __name__ == "__main__":
    # create parser
    parser = argparse.ArgumentParser()
    # define args
    parser.add_argument("-b", "--bed_file", type=Path)
    parser.add_argument("-g", "--genome_file", type=Path)
    parser.add_argument("-o", "--output_file", type=Path)
    parser.add_argument("-s", "--surrounding_bases", type=int, default=1)
    parser.add_argument("--motif_col_position", choices=(0, 1), default=0)
    parser.add_argument("--keep_temp_files", action="store_true")
    parser.add_argument("--skip_header", action="store_true")
    # parse args
    args = parser.parse_args()
    # run script with parsed args
    single_logo(
        args.bed_file,
        args.genome_file,
        args.output_file,
        args.surrounding_bases,
        args.motif_col_position,
        args.keep_temp_files,
        args.skip_header,
    )

