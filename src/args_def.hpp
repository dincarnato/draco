#pragma once

#include "args_impl.hpp"
#include <string>

namespace args {

static constexpr auto opts = args::Opts(
    "\n DRACO (v1.0)\n",
    args::Group(
        "",
        ARG(std::string, mm_filename)
            .parameter_name("mm")
            .description("Input mutation map (MM) file"),
        ARG(std::string, output_filename)
            .parameter_name("output")
            .description("Output JSON file")
            .DEFAULT_VALUE("draco_deconvoluted.json"),
        ARG(unsigned, n_processors)
            .parameter_name("processors")
            .description("Number of processors to use [Note: when set to 0, "
                         "all available processors will be used. Only the analysis of multiple "
                         "transcripts can be parallelized; if analyzing a single transcript, "
                         "setting this parameter to a value > 1 will not speed up the execution]")
            .DEFAULT_VALUE(0),
        ARG(std::string, whitelist)
            .optional()
            .parameter_name("whitelist")
            .description("A whitelist file, containing the IDs of the transcripts "
                         "to be analyzed, one per row"),
        ARG(bool, shape)
            .description(
                "Enables spectral analysis on all four bases (default is only A/C bases) "
                "[Note: this feature is highly experimental]")
            .DEFAULT_VALUE(false)),

    args::Group("Mutation filtering",
        ARG(unsigned, minimum_base_coverage)
            .parameter_name("minBaseCoverage")
            .description("Minimum coverage per base")
            .DEFAULT_VALUE(100u),
        ARG(unsigned, minimum_modifications_per_base)
            .parameter_name("minBaseMutations")
            .description("Minimum mutations per base")
            .DEFAULT_VALUE(2u),
        ARG(float, minimum_modifications_per_base_fraction)
            .parameter_name("minMutationFreq")
            .description("Bases below this mutation frequency will be discarded as noise")
            .DEFAULT_VALUE(0.005f),
        ARG(unsigned, minimum_modifications_per_read)
            .parameter_name("minReadMutations")
            .description("Reads with less than these mutations will be discarded as non-informative")
            .DEFAULT_VALUE(2u)),

    args::Group(
        "Spectral deconvolution",
        ARG(unsigned, max_clusters)
            .parameter_name("maxClusters")
            .description("Maximum allowed number of clusters/conformations")
            .DEFAULT_VALUE(std::numeric_limits<unsigned>::max()),
        ARG(unsigned, min_filtered_reads)
            .parameter_name("minFilteredReads")
            .description("Minimum number of reads (post-filtering) to perform spectral analysis")
            .DEFAULT_VALUE(5),
        ARG(unsigned, min_permutations)
            .parameter_name("minPermutations")
            .description("Minimum number of permutations performed to build the null model")
            .DEFAULT_VALUE(8),
        ARG(unsigned, max_permutations)
            .parameter_name("maxPermutations")
            .description("Maximum number of permutations performed to build the null model")
            .DEFAULT_VALUE(400),
        ARG(double, first_eigengap_threshold)
            .parameter_name("firstEigengapThresh")
            .description("Threshold to consider the first eigengap [Note: when this threshold "
                         "is not met, 0 clusters are reported]")
            .DEFAULT_VALUE(0.90),
        ARG(double, min_eigengap_threshold)
            .parameter_name("eigengapCumRelThresh")
            .description("Minimum relative difference between the eigengap and the null model, "
                         "as a fraction of the cumulative difference between the previous eigengaps "
                         "and their respective null models [Note: this does not apply to the first eigengap]")
            .DEFAULT_VALUE(0.10),
        ARG(double, eigengap_diff_absolute_threshold)
            .parameter_name("eigengapAbsThresh")
            .description("Minimum absolute difference between the eigengap and the null model")
            .DEFAULT_VALUE(0.03),
        ARG(double, alpha_value)
            .parameter_name("alpha")
            .description("Below this p-value, the null hypothesis is rejected "
                         "and the eigengap is marked as informative")
            .DEFAULT_VALUE(0.01),
        ARG(double, beta_value)
            .parameter_name("beta")
            .description("Above this p-value, the alternative hypothesis is rejected "
                         "and the eigengap is marked as non-informative "
                         "[Note: this threshold does not apply to the first eigengap]")
            .DEFAULT_VALUE(0.2),
        ARG(double, first_eigengap_beta_value)
            .parameter_name("firstEigengapBeta")
            .description("Beta p-value threshold for the first eigengap")
            .DEFAULT_VALUE(0.4),
        ARG(double, min_null_stddev)
            .parameter_name("minNullStdev")
            .description("Minimum standard deviation for the null model "
                         "[Note: when this threshold is not met, \"--extraPermutations\" additional permutations "
                         "will be performed]")
            .DEFAULT_VALUE(0.025),
        ARG(unsigned, alternative_check_permutations)
            .parameter_name("extraPermutations")
            .description("Additional permutations to perform when the standard deviation of the null "
                         " model is < \"--minNullStdev\"")
            .DEFAULT_VALUE(50),
        ARG(unsigned, min_bases_size)
            .parameter_name("minWinBases")
            .description("Minimum number of bases in window (post-filtering) to perform spectral analysis")
            .DEFAULT_VALUE(10),
        ARG(unsigned char, extended_search_eigengaps)
            .parameter_name("lookaheadEigengaps")
            .description("Number of eigengaps to look ahead after a non-informative eigengap is encountered")
            .DEFAULT_VALUE(3),
        ARG(bool, create_eigengaps_plots)
            .parameter_name("saveEigengapData")
            .description("Saves eigengap data for plotting")
            .DEFAULT_VALUE(false),
        ARG(std::string, eigengaps_plots_root_dir)
            .parameter_name("eigengapDataOut")
            .description("Eigengap data output folder")
            .DEFAULT_VALUE("./eigengap_data")),

    args::Group(
        "Graph-Cut",
        ARG(double, minimum_cluster_fraction)
            .parameter_name("minClusterFraction")
            .description("Minimum fraction of reads assigned to each cluster/conformation "
                         "[Note: if this threshold is not met, the number of clusters is automatically decreased]")
            .DEFAULT_VALUE(0.05)),

    args::Group(
        "Windowed analysis",
        ARG(double, window_size_fraction)
            .parameter_name("winLenFraction")
            .description("Length of the window as fraction of the median read length [Note: this parameter "
                         "and \"--absWinSize\" are mutually exclusive")
            .DEFAULT_VALUE(0.9),
        ARG(unsigned, window_size)
            .parameter_name("absWinLen")
            .description("Absolute length of the window [Note: this parameter and \"--winSizeFraction\" "
                         "are mutually exclusive")
            .DEFAULT_VALUE(0u),
        ARG(double, window_shift_fraction)
            .parameter_name("winOffsetFraction")
            .description("Slide offset as fraction of the size of the window [Note: this parameter and \"--absWinOffset\" "
                         "are mutually exclusive")
            .DEFAULT_VALUE(0.05),
        ARG(unsigned, window_shift)
            .parameter_name("absWinOffset")
            .description("Absolute slide offset [Note: this parameter and \"--winOffsetFraction\" "
                         "are mutually exclusive")
            .DEFAULT_VALUE(0u),
        ARG(unsigned, max_collapsing_windows)
            .parameter_name("maxIgnoreWins")
            .description("Maximum number of internal windows with a different number of "
                         "clusters to ignore when merging two external sets of windows")
            .DEFAULT_VALUE(0),
        ARG(unsigned, min_surrounding_windows_size)
            .parameter_name("minExtWins")
            .description("Minimum number of external windows, having the same number of clusters, "
                         "needed to trigger merging")
            .DEFAULT_VALUE(6),
        ARG(bool, set_uninformative_clusters_to_surrounding)
            .parameter_name("nonInformativeToSurround")
            .description("Non-informative windows (windows with 0 detected clusters) are set to the same "
                         "number of clusters of surrounding windows")
            .DEFAULT_VALUE(false),
        ARG(bool, set_all_uninformative_to_one)
            .parameter_name("allNonInformativeToOne")
            .description("If all windows in the transcript are non-informative (0 clusters), the number of clusters is "
                         "forced to 1")
            .DEFAULT_VALUE(false),
        ARG(bool, report_uninformative)
            .parameter_name("reportNonInformative")
            .description("Reports also non-informative windows in the output JSON file")
            .DEFAULT_VALUE(false)));
} // namespace args
