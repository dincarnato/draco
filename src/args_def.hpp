#pragma once

#include "args_impl.hpp"
#include <string>

namespace args {

static constexpr auto opts = args::Opts(
    "Analyze an MM file, evaluating the structural ensemble for each "
    "transcript",
    args::Group(
        "",
        ARG(std::string, mm_filename)
            .parameter_name("mm")
            .description("Input MM filename"),
        ARG(std::string, output_filename)
            .parameter_name("output")
            .description("Output filename")
            .DEFAULT_VALUE("out.json"),
        ARG(unsigned, n_processors)
            .description("Number of processors to use. When set to 0, the "
                         "maximum number of processors will be used.")
            .DEFAULT_VALUE(0),
        ARG(std::string, whitelist)
            .optional()
            .parameter_name("whitelist")
            .description("A filename with a whitelist of genes. Each gene must "
                         "be in a separate line."),
        ARG(bool, shape)
            .DEFAULT_VALUE(false)
            .description(
                "A SHAPE reagent is used in the analysis instead of DMS "
                "(default). This causes all the bases to be considered as "
                "possibly modified instead of A and C only.")),
    args::Group("Modifications matrix filtering",
                ARG(unsigned, minimum_base_coverage).DEFAULT_VALUE(100u),
                ARG(unsigned, minimum_modifications_per_base).DEFAULT_VALUE(2u),
                ARG(unsigned, minimum_modifications_per_read).DEFAULT_VALUE(2u),
                ARG(float, minimum_modifications_per_base_fraction)
                    .DEFAULT_VALUE(0.005f)),

    args::Group(
        "Number of clusters detection",
        ARG(unsigned, min_filtered_reads).DEFAULT_VALUE(5),
        ARG(unsigned, max_permutations).DEFAULT_VALUE(400),
        ARG(unsigned, min_permutations).DEFAULT_VALUE(8),
        ARG(double, first_eigengap_threshold).DEFAULT_VALUE(0.90),
        ARG(double, min_eigengap_threshold).DEFAULT_VALUE(0.10),
        ARG(double, eigengap_diff_absolute_threshold).DEFAULT_VALUE(0.03),
        ARG(double, alpha_value).DEFAULT_VALUE(0.01),
        ARG(double, beta_value).DEFAULT_VALUE(0.2),
        ARG(double, first_eigengap_beta_value).DEFAULT_VALUE(0.4),
        ARG(unsigned, max_clusters)
            .DEFAULT_VALUE(std::numeric_limits<unsigned>::max()),
        ARG(unsigned, alternative_check_permutations).DEFAULT_VALUE(50),
        ARG(double, min_null_stddev).DEFAULT_VALUE(0.025),
        ARG(unsigned, min_bases_size).DEFAULT_VALUE(10),
        ARG(unsigned char, extended_search_eigengaps).DEFAULT_VALUE(3),
        ARG(bool, create_eigengaps_plots).DEFAULT_VALUE(false),
        ARG(std::string, eigengaps_plots_root_dir)
            .DEFAULT_VALUE("./eigengaps_plots")),

    args::Group(
        "Results refinement",
        ARG(double, minimum_cluster_fraction)
            .DEFAULT_VALUE(0.05)
            .parameter_name("min-cluster-fraction")
            .description("The minimum fraction of reads assigned to each "
                         "cluster. When a fraction of reads assigned to a "
                         "cluster is below this threshold, the number of "
                         "clusters is automatically decreased.")),

    args::Group(
        "Windows collapsing",
        ARG(double, window_size_fraction)
            .DEFAULT_VALUE(0.9)
            .description(
                "The size of the window as fraction of the median size of "
                "the reads"),
        ARG(unsigned, window_size)
            .DEFAULT_VALUE(0u)
            .description(
                "The size of the window. If this parameter is specified, the "
                "'window_size_fraction' parameter is ignored"),
        ARG(double, window_shift_fraction)
            .DEFAULT_VALUE(0.05)
            .description(
                "The shift of the window as fraction of the window size"),
        ARG(unsigned, window_shift)
            .DEFAULT_VALUE(0u)
            .description(
                "The shift of the window. If this parameter is specified, the "
                "'window_shift_fraction' parameter is ignored"),
        ARG(unsigned, max_collapsing_windows).DEFAULT_VALUE(0),
        ARG(unsigned, min_surrounding_windows_size).DEFAULT_VALUE(6),
        ARG(bool, set_uninformative_clusters_to_surrounding)
            .DEFAULT_VALUE(false),
        ARG(bool, set_all_uninformative_to_one)
            .DEFAULT_VALUE(false)
            .description(
                "In case all the windows are uninformative, the number of "
                "clusters for each window is set to 1"),
        ARG(bool, report_uninformative)
            .DEFAULT_VALUE(false)
            .description("Report the uninformative windows into the report")));
} // namespace args
