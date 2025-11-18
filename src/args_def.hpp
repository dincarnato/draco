#pragma once

#include "args_impl.hpp"

#include <string>

namespace args {

static constexpr auto opts = args::Opts(
    args::Group(
        "",
        ARG(std::string, mm_filenames)
            .parameter_name("mm")
            .multiplicity<multiplicity::Many>()
            .description("Input mutation map (MM) file [Note: to pass "
                         "multiple replicates, simply call the parameter "
                         "multiple times (e.g., --mm rep1.mm --mm rep2.mm)]"),
        ARG(std::string, output_filename)
            .parameter_name("output")
            .description("Output JSON file")
            .DEFAULT_VALUE("draco_deconvolved.json"),
        ARG(unsigned, n_processors)
            .parameter_name("processors")
            .description("Number of processors to use [Note: when set to 0, "
                         "all available processors will be used]")
            .DEFAULT_VALUE(0),
        ARG(std::string, whitelist)
            .optional()
            .parameter_name("whitelist")
            .description(
                "A whitelist file, containing the IDs of the transcripts "
                "to be analyzed, one per row"),
        ARG(bool, shape)
            .description(
                "Enables spectral analysis on all four bases (Default: "
                "only A/C bases) [Note: this feature is highly experimental]")
            .DEFAULT_VALUE(false),
        ARG(std::string, output_raw_n_clusters)
            .optional()
            .parameter_name("outputRawNClusters")
            .description(
                "Outputs a BED-like file containing, for each transcript, the "
                "raw number of clusters detected for each window"),
        ARG(std::string, log_level)
            .optional()
            .parameter_name("log-level")
            .description("Specifies the log level (allowed values: trace, "
                         "debug, info, warn, error)")),

    args::Group(
        "Mutation filtering",
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
            .description("Bases below this mutation frequency will be "
                         "discarded as noise")
            .DEFAULT_VALUE(0.005f),
        ARG(float, maximum_modifications_per_base_fraction)
            .parameter_name("maxMutationFreq")
            .description("Bases above or equal this mutation frequency will be "
                         "discarded as noise")
            .DEFAULT_VALUE(0.35f),
        ARG(unsigned, minimum_modifications_per_read)
            .parameter_name("minReadMutations")
            .description(
                "Reads with fewer than this number of mutations will be "
                "discarded as non-informative")
            .DEFAULT_VALUE(2u)),

    args::Group(
        "Spectral deconvolution",
        ARG(unsigned, max_clusters)
            .parameter_name("maxClusters")
            .description("Maximum allowed number of clusters/conformations")
            .DEFAULT_VALUE(std::numeric_limits<unsigned>::max()),
        ARG(unsigned, min_filtered_reads)
            .parameter_name("minFilteredReads")
            .description("Minimum number of reads (post-filtering) to perform "
                         "spectral analysis")
            .DEFAULT_VALUE(100),
        ARG(unsigned, min_permutations)
            .parameter_name("minPermutations")
            .description("Minimum number of permutations performed to build "
                         "the null model")
            .DEFAULT_VALUE(10),
        ARG(unsigned, max_permutations)
            .parameter_name("maxPermutations")
            .description("Maximum number of permutations performed to build "
                         "the null model")
            .DEFAULT_VALUE(100),
        ARG(bool, ignore_first_eigengap)
            .parameter_name("ignoreFirstEigengap")
            .description("The first eigengap is ignored")
            .DEFAULT_VALUE(false),
        ARG(double, min_eigengap_threshold)
            .parameter_name("eigengapCumRelThresh")
            .description("Minimum relative difference between the eigengap and "
                         "the null model, "
                         "as a fraction of the cumulative difference between "
                         "the previous eigengaps "
                         "and their respective null models")
            .DEFAULT_VALUE(0.10),
        ARG(double, eigengap_diff_absolute_threshold)
            .parameter_name("eigengapAbsThresh")
            .description("Minimum absolute difference between the eigengap and "
                         "its null model")
            .DEFAULT_VALUE(0.03),
        ARG(double, alpha_value)
            .parameter_name("alpha")
            .description("Below this p-value, the null hypothesis is rejected "
                         "and the eigengap is marked as informative")
            .DEFAULT_VALUE(0.01),
        ARG(double, beta_value)
            .parameter_name("beta")
            .description(
                "Above this p-value, the alternative hypothesis is rejected "
                "and the eigengap is marked as non-informative")
            .DEFAULT_VALUE(0.2),
        ARG(double, min_null_stddev)
            .parameter_name("minNullStdev")
            .description("Minimum standard deviation for the null model "
                         "[Note: when this threshold is not met, "
                         "\"--extraPermutations\" additional permutations "
                         "will be performed]")
            .DEFAULT_VALUE(0.025),
        ARG(unsigned, alternative_check_permutations)
            .parameter_name("extraPermutations")
            .description("Additional permutations to perform when the standard "
                         "deviation of the null "
                         " model is < \"--minNullStdev\"")
            .DEFAULT_VALUE(50),
        ARG(unsigned, min_bases_size)
            .parameter_name("minWinBases")
            .description("Minimum number of bases in window (post-filtering) "
                         "to perform spectral analysis")
            .DEFAULT_VALUE(10),
        ARG(unsigned char, extended_search_eigengaps)
            .parameter_name("lookaheadEigengaps")
            .description("Number of eigengaps to look ahead after a "
                         "non-informative eigengap is encountered")
            .DEFAULT_VALUE(0),
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
            .description("Minimum fraction of reads assigned to each "
                         "cluster/conformation "
                         "[Note: if this threshold is not met, the number of "
                         "clusters is automatically decreased]")
            .DEFAULT_VALUE(0.005),
        ARG(std::uint16_t, soft_clustering_kmeans_iterations)
            .parameter_name("softClusteringKmeansIters")
            .description("Number of iterations of kmeans performed on "
                         "graph-cut when more than 2 clusters are found")
            .DEFAULT_VALUE(20)),

    args::Group(
        "Windowed analysis",
        ARG(double, window_size)
            .parameter_name("winLen")
            .description(
                "Length of the window. If this value is comprised between 0 "
                "and 1, it "
                "is interpreted as a fraction of the median read length. If "
                "this value "
                "is > 1, it is interpreted as the absolute length of the "
                "window. [Note: this "
                "parameter and \"--winLenFracRnaLen\"  are mutually exclusive]")
            .DEFAULT_VALUE(100),
        ARG(double, window_size_fraction_transcript_size)
            .parameter_name("winLenFracRnaLen")
            .description("Length of the window as fraction of length of the "
                         "RNA [Note: this parameter and \"--winLen\" are "
                         "mutually exclusive]")
            .DEFAULT_VALUE(0),
        ARG(double, window_shift)
            .parameter_name("winOffset")
            .description("Window sliding offset. If this value is comprised "
                         "between 0 and 1, it is "
                         "interpreted as a fraction of the window's length. If "
                         "this value is >= 1, "
                         "it is interpreted as the number of bases to slide "
                         "the window by.")
            .DEFAULT_VALUE(0.01),
        ARG(unsigned, max_collapsing_windows)
            .parameter_name("maxIgnoreWins")
            .description(
                "Maximum number of internal windows with a different number of "
                "clusters to ignore when merging two external sets of windows")
            .DEFAULT_VALUE(0),
        ARG(unsigned, min_surrounding_windows_size)
            .parameter_name("minExtWins")
            .description("Minimum number of external windows, having the same "
                         "number of clusters, needed to trigger merging")
            .DEFAULT_VALUE(6),
        ARG(bool, set_uninformative_clusters_to_surrounding)
            .parameter_name("nonInformativeToSurround")
            .description("Non-informative windows (windows with 0 detected "
                         "clusters) are set to the same "
                         "number of clusters of surrounding windows")
            .DEFAULT_VALUE(false),
        ARG(bool, set_all_uninformative_to_one)
            .parameter_name("allNonInformativeToOne")
            .description("If all windows in the transcript are non-informative "
                         "(0 clusters), the number of clusters is forced to 1")
            .DEFAULT_VALUE(false),
        ARG(bool, report_uninformative)
            .parameter_name("reportNonInformative")
            .description("Reports also non-informative windows (0 clusters) in "
                         "the output JSON file")
            .DEFAULT_VALUE(false),
        ARG(std::string, assignments_dump_directory)
            .optional()
            .parameter_name("assignmentsDumpDir")
            .description("When specified, a folder is generated containing an "
                         "MM file for "
                         "each set of merged windows, containing a dump of the "
                         "reads assigned "
                         "to each conformation"),
        ARG(bool, skip_ambiguous_assignments)
            .parameter_name("skipAmbiguousAssignments")
            .description(
                "When specified, if the best assignment score for a read is "
                "equal across multiple clusters, the read is discarded")
            .DEFAULT_VALUE(false),
        ARG(float, min_windows_overlap)
            .parameter_name("minWindowsOverlap")
            .description(
                "Minimum overlap between two non-contiguous "
                "windows with the same number of conformations in "
                "order to merge them in the same region. The value must be "
                "comprised between 0 and 1, and it is interpreted as a "
                "fraction "
                "of the window size")
            .DEFAULT_VALUE(1.f)));
} // namespace args
