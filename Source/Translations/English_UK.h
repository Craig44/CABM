//============================================================================
// Name        : English.h
// Author      : S.Rasmussen
// Date        : 18/07/2008
// Copyright   : Copyright NIWA Science �2008 - www.niwa.co.nz
// Description :
// $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
//============================================================================
#ifndef ENGLISH_UK_H_
#define ENGLISH_UK_H_

// WARNING
// TEXT STRING DEFINED AS PARAM_X MUST BE LOWERCASE ONLY

//**********************************************************************
// BASE CONFIGURATION
//
//**********************************************************************
// Configuration File Syntax Strings
#define CONFIG_ARRAY_END                    "]"
#define CONFIG_ARRAY_START                  "["
#define CONFIG_CATEGORY_SEPARATOR           "."
#define CONFIG_END_REPORT                   "*end"
#define CONFIG_FALSE                        "false"
#define CONFIG_FALSE_SHORT                  "f"
#define CONFIG_INCLUDE                      "!include"
#define CONFIG_JOIN_OPERATOR                "+"
#define CONFIG_LIST_OPERATOR                ","
#define CONFIG_MULTI_COMMENT_START          "/*"
#define CONFIG_MULTI_COMMENT_END            "*/"
#define CONFIG_RATIO_SEPARATOR              ":"
#define CONFIG_RANGE_OPERATOR               ":"
#define CONFIG_SECTION_SYMBOL               "@"
#define CONFIG_SEPERATOR_ESTIMATE_VALUES    ","
#define CONFIG_SINGLE_COMMENT               "#"
#define CONFIG_TRUE                         "true"
#define CONFIG_TRUE_SHORT                   "t"
#define CONFIG_WILDCARD_MULTIPLE            "*"
#define CONFIG_WILDCARD_SINGLE              "?"

//**********************************************************************
// REPORT
//
//**********************************************************************
// Report output Syntax Strings
#define REPORT_END                           "*end"
#define REPORT_R_COMPLETE_VECTOR             "{c}"    // not sure if this is used
#define REPORT_R_DATAFRAME                   "{d}"    // R library will add a header to the data frame, so if you don't add one to the report R will put row 1 as the header
#define REPORT_R_MATRIX                      "{m}"    // No header's
#define REPORT_R_NAMED_COMPLETE_VECTOR       "{C}"    // not sure if this is used
#define REPORT_R_LIST                        "{L}"
#define REPORT_R_LIST_END                    "end {L}"
#define REPORT_R_LIST_ELEMENT_SEPERATOR      ":"
#define REPORT_R_VECTOR                      "{v}"
#define REPORT_R_STRING_VECTOR               "{s}"    // This will not split a vector: used to read in warnings, see Logging.cpp warnings section

//**********************************************************************
// PARAMETERS
//
//**********************************************************************
#define PARAM_A                                   "a"
#define PARAM_A_LAYER_LABEL                       "a_layer_label"
#define PARAM_A50                                 "a50"
#define PARAM_ABUNDANCE                           "abundance"
#define PARAM_ABUNDANCE_DENSITY                   "abundance_density"
#define PARAM_ACTIVE                              "active"
#define PARAM_ADDITIVE                            "additive"
#define PARAM_ADDITIVE_NORMAL                     "additive_normal"
#define PARAM_ADDRESSABLE                         "addressable"
#define PARAM_ADJACENT_CELL_MOVEMENT              "adjacent_cell"
#define PARAM_AGE                                 "age"
#define PARAM_AGE_WEIGHT_LABELS                   "age_weight_labels"
#define PARAM_AGES                                "ages"
#define PARAM_AGE_FREQUENCY_BY_CELL               "age_frequency_by_cell"
#define PARAM_AGE_INDEX                           "age_index"
#define PARAM_AGE_LENGTH                          "age_length"
#define PARAM_AGE_LENGTHS                         "age_lengths"
#define PARAM_AGE_LENGTH_KEY                      "age_length_key"
#define PARAM_AGE_LENGTH_KEYS                     "age_length_keys"
#define PARAM_AGE_PLUS                            "age_plus"
#define PARAM_AGE_V                               "age_v"
#define PARAM_AGE_WEIGHT                          "age_weight"
#define PARAM_AGE_WEIGHT_LABEL                    "age_weight_label"
#define PARAM_AGEING                              "ageing"
#define PARAM_AGEING_LABEL                        "ageing_label"
#define PARAM_AGEING_ERROR                        "ageing_error"
#define PARAM_AGEING_ERRORS                       "ageing_errors"
#define PARAM_AGENTS							                "agents"
#define PARAM_AGENT                               "agent"
#define PARAM_ALL_VALUES_BOUNDED                  "all_values_bounded"
#define PARAM_ALL_VALUES                          "all_values"
#define PARAM_ALPHA                               "alpha"
#define PARAM_ANNUAL_MORTALITY_RATE               "annual_mortality_rate"
#define PARAM_ANNUAL_SHIFT                        "annual_shift"
#define PARAM_APPEND                              "append"
#define PARAM_AREA                                "area"
#define PARAM_AREAS                               "areas"
#define PARAM_ARMA                                "arma"
#define PARAM_ASSERT                              "assert"
#define PARAM_ATO95                               "ato95"
#define PARAM_AVERAGE                             "average"
#define PARAM_AVERAGE_DIFFERENCE                  "average_difference"
#define PARAM_AVERAGE_LOWER_BOUND                 "average_lower_bound"
#define PARAM_AVERAGE_UPPER_BOUND                 "average_upper_bound"
#define PARAM_B                                   "b"
#define PARAM_B_LAYER_LABEL                       "b_layer_label"
#define PARAM_BASE_UNTIS                          "base_weight_units"
#define PARAM_BASE_LAYER_LABEL                    "base_layer_label"
#define PARAM_B0                                  "b0"
#define PARAM_B0_VALUE                            "b0_value"
#define PARAM_B0_PHASE                            "b0_intialisation_phase"
#define PARAM_BASIC                               "basic"
#define PARAM_BETA                                "beta"
#define PARAM_BEVERTON_HOLT                       "beverton_holt"
#define PARAM_BH_RECRUITMENT                      "bh_recruitment"
#define PARAM_BINOMIAL                            "binomial"
#define PARAM_BINOMIAL_APPROX                     "binomial_approx"
#define PARAM_BIN_LABELS                          "bin_labels"
#define PARAM_BIOMASS                             "biomass"
#define PARAM_BIOMASS_DENSITY                     "biomass_density"
#define PARAM_BIOMASS_EVENT_MORTALITY             "biomass_event_mortality"
#define PARAM_BIOMASS_LAYER_LABEL                 "biomass_layer_label"
#define PARAM_BY_LENGTH                           "by_length"
#define PARAM_C                                   "c"
#define PARAM_CATCHABILITY                        "catchability"
#define PARAM_CATCH                               "catch"
#define PARAM_CATCH_LAYERS                        "catch_layers"
#define PARAM_CATCHES                             "catches"
#define PARAM_CATEGORICAL                         "categorical"
#define PARAM_CATEGORICAL_MONOTONIC               "monotonic_categorical"
#define PARAM_CATEGORIES                          "categories"
#define PARAM_CATEGORY                            "category"
#define PARAM_CELL_LENGTH                         "cell_length"
#define PARAM_CELLS                               "cells"
#define PARAM_CLASS_MINIMUMS                      "class_minimums"
#define PARAM_COLUMN                              "column"
#define PARAM_COLUMN_INDEX                        "column_index"
#define PARAM_CONFIG_FILE                         "config_file"
#define PARAM_CONSTANT                            "constant"
#define PARAM_CONSTANT_LOSS_RATE                  "constant_loss_rate"
#define PARAM_CONSTANT_RATE                       "constant_rate"
#define PARAM_CONSTANT_RECRUITMENT                "constant_recruitment"
#define PARAM_CONSUMPTION_RATE                    "consumption_rate"
#define PARAM_CONVERGENCE_YEARS                   "convergence_years"
#define PARAM_CROSSOVER_PROBABILITY               "crossover_probability"
#define PARAM_CURRENT_YEAR                        "current_year"
#define PARAM_CV                                  "cv"
#define PARAM_D                                   "d"
#define PARAM_DATA                                "data"
#define PARAM_DATA_WEIGHT_TYPE                    "data_weight_type"
#define PARAM_DATA_WEIGHT_VALUE                   "data_weight_value"
#define PARAM_DE_SOLVER                           "de_solver"
#define PARAM_DEBUG                               "debug"
#define PARAM_DELTA                               "delta"
#define PARAM_DEFAULT_LAYER                       "default_layer"
#define PARAM_DERIVED                             "derived"
#define PARAM_DERIVED_QUANTITIES                  "derived_quantities"
#define PARAM_DERIVED_QUANTITY                    "derived_quantity"
#define PARAM_DETECTION_PARAMETER                 "detection"
#define PARAM_DEVIATION_VALUES                    "deviation_values"
#define PARAM_DEVIATION_YEARS                     "deviation_years"
#define PARAM_DF                                  "df"
#define PARAM_DIFFERENCE                          "difference"
#define PARAM_DISPERSION                          "dispersion"
#define PARAM_DISTANCE                            "distance"
#define PARAM_DISCARD_MORTALITY_RATE              "discard_mortality_rate"
#define PARAM_DISTRIBUTION                        "distribution"
#define PARAM_DIRICHLET					            		  "dirichlet"
#define PARAM_DOUBLE                              "numeric"
#define PARAM_DOUBLE_EXPONENTIAL                  "double_exponential"
#define PARAM_DOUBLE_NORMAL                       "double_normal"
#define PARAM_E                                   "e"
#define PARAM_ELECTIVITIES                        "electivities"
#define PARAM_ELEMENT_DIFFERENCE                  "element_difference"
#define PARAM_EMPIRICAL_SAMPLING                  "empirical_sampling"
#define PARAM_ENABLED_ESTIMATES                   "enabled_estimates"
#define PARAM_END_TABLE                           "end_table"
#define PARAM_EQUATION                            "equation"
#define PARAM_EQUILIBRIUM_METHOD                  "equilibrium_method"
#define PARAM_ERROR_VALUE                         "error_value"
#define PARAM_ERROR_VALUES                        "error_values"
#define PARAM_ERROR_VALUE_MULTIPLIER              "error_value_multiplier"
#define PARAM_EVENT                               "event"
#define PARAM_EVENT_MORTALITY                     "event_mortality"
#define PARAM_EXCLUDE_PROCESSES                   "exclude_processes"
#define PARAM_EXPECTED_VALUE                      "expected_value"
#define PARAM_EXPONENTIAL                         "exponential"
#define PARAM_EXOGENOUS_VARIABLE                  "exogeneous_variable"
#define PARAM_EXOGENOUS                           "exogeneous"
#define PARAM_FILE_NAME                           "file_name"
#define PARAM_FIRST_YEAR                          "first_year"
#define PARAM_FINAL_YEAR                          "final_year"
#define PARAM_FINE                                "fine"
#define PARAM_FINEST                              "finest"
#define PARAM_FISHERY                             "fishery"
#define PARAM_FISHERY_TIME_STEPS                  "fishery_time_step"
#define PARAM_FISHERIES                           "fisheries"
#define PARAM_FORCE_ESTIMABLE_VALUES_FILE         "force_estimable_values_file"
#define PARAM_FORMAT                              "format"
#define PARAM_FRANCIS                             "francis"
#define PARAM_FREE                                "free"
#define PARAM_FROM                                "from"
#define PARAM_FUNCTION                            "function"
#define PARAM_GRAMS                               "grams"
#define PARAM_GROWTH                              "growth"
#define PARAM_GROWTH_BASED                        "growth_based"
#define PARAM_GROWTH_BASIC                        "growth_basic"
#define PARAM_GROWTH_PROPORTIONS                  "growth_proportions"
#define	PARAM_GROWTH_PROCESS_LABEL				        "growth_process_label"
#define PARAM_GROWTH_TIME_STEPS                   "growth_time_steps"
#define PARAM_GROWTH_VON_BERTALANFFY              "growth_von_bertalanffy"
#define PARAM_GROWTH_VON_BERTALANFFY_WITH_BASIC   "growth_von_bertalanffy_with_basic_weight"
#define PARAM_H                                   "h"
#define PARAM_HEADER                              "header"
#define PARAM_HEIGHT                              "height"
#define PARAM_HYBRID                              "hybrid"
#define PARAM_INCREASING                          "increasing"
#define PARAM_INDEX                               "index"
#define PARAM_INITIAL_YEAR                        "initial_year"
#define PARAM_INITIALISATION                      "initialisation"
#define PARAM_INITIALISATION_TIME_STEPS           "initialisation_time_steps"
#define PARAM_INITIALISATION_PARTITION            "initialisation_partition"
#define PARAM_INITIALISATION_PHASE_LABELS          "initialisation_phase_labels"
#define PARAM_INITIALISATION_PHASE                 "initialisation_phase"
#define PARAM_INITIALISATION_PHASES               "initialisation_phases"
#define PARAM_INITIALISATION_SEED_Z               "initialisation_z_seed"
#define PARAM_INSERT_PROCESSES                    "insert_processes"
#define PARAM_INTEGER                             "integer"
#define PARAM_INTERVALS                           "intervals"
#define PARAM_IS_ABUNDANCE                        "is_abundance"
#define PARAM_INCREMENTAL_SUFFIX                  "incremental_suffix"
#define PARAM_INTEGER                             "integer"
#define PARAM_INTERNAL_GAPS                       "internal_gaps"
#define PARAM_INTERPOLATE                         "interpolate"
#define PARAM_INTERCEPT                           "intercept"
#define PARAM_INVERSE                             "inverse"
#define PARAM_INVERSE_LOGISTIC                    "inverse_logistic"
#define PARAM_ITERATIVE                           "iterative"
#define PARAM_JACOBIAN                            "jacobian"
#define PARAM_K                                   "k"
#define PARAM_K_LAYER_LABEL					              "k_layer_label"
#define PARAM_KGS                                 "kgs"
#define PARAM_KEEP                                "keep"
#define PARAM_KNIFE_EDGE                          "knife_edge"
#define PARAM_L                                   "l"
#define PARAM_LABEL                               "label"
#define PARAM_LAMBDA                              "lambda"
#define PARAM_LAST_YEAR_WITH_NO_BIAS              "last_year_with_no_bias"
#define PARAM_LAST_YEAR_WITH_BIAS                 "last_year_with_bias"
#define PARAM_LATITUDE_LAYER_LABEL                "latitude_layer_label"
#define PARAM_LAYER                               "layer"
#define PARAM_LAYER_LABELS                        "layer_labels"
#define PARAM_LAYERS                              "layers"
#define PARAM_LAYER_DERIVED_WORLD_VIEW            "layer_derived_view"
#define PARAM_LAYER_HEIGHT                        "layer_height"
#define PARAM_LAYER_LABEL                         "layer_label"
#define PARAM_LAYER_WIDTH                         "layer_width"
#define PARAM_LAYERS                              "layers"
#define PARAM_LAYER_OF_CELLS                      "cell_layer"
#define PARAM_LENGTH                              "length"
#define PARAM_LENGTH_BASED                        "length_based"
#define PARAM_LENGTH_BINS                         "length_bins"
#define PARAM_LENGTH_PLUS                         "length_plus"
#define PARAM_LENGTH_V                            "length_v"
#define PARAM_LENGTH_VALUES                       "length_values"
#define PARAM_LENGTH_WEIGHT                       "length_weight"
#define PARAM_LENGTH_WEIGHTS                      "length_weights"
#define PARAM_LIKELIHOOD                          "likelihood"
#define PARAM_LIKELIHOOD_MULTIPLIER               "likelihood_multiplier"
#define PARAM_LINF                                "linf"
#define PARAM_LINF_LAYER_LABEL                    "linf_layer_label"
#define PARAM_LINEAR_INTERPOLATION                "linear_interpolation"
#define PARAM_LINEAR                              "linear"
#define PARAM_LOCAL_BH_RECRUITMENT                "local_bh_recruitment"
#define PARAM_LOG                                 "log"
#define PARAM_LOG_LEVEL                           "log_level"
#define PARAM_LOG_ODDS                            "log_odds"
#define PARAM_LOG_SCALE                           "log_scale"
#define PARAM_LOG_SUM                             "log_sum"
#define PARAM_LOGISTIC                            "logistic"
#define PARAM_LOGISTIC_NORMAL                     "logistic_normal"
#define PARAM_LOGISTIC_PRODUCING                  "logistic_producing"
#define PARAM_LOGNORMAL                           "lognormal"
#define PARAM_LOGNORMAL_EMPIRICAL                 "lognormal_empirical"
#define PARAM_LOGNORMAL_TRANSFORMATION            "lognormal_transfomration"
#define PARAM_LOGNORMAL_WITH_Q                    "lognormal_with_q"
#define PARAM_LONGITUDE_LAYER_LABEL               "longitude_layer_label"
#define PARAM_LOSS_RATE                           "loss_rate"
#define PARAM_LOSS_RATE_SELECTIVITIES             "loss_rate_selectivities"
#define PARAM_LOWER_BOUND                         "lower_bound"
#define PARAM_M                                   "m"
#define PARAM_M_MULTIPLIER_LAYER_LABEL            "m_multiplier_layer_label"
#define PARAM_MARKOVIAN                           "markovian"
#define PARAM_MATURE                              "mature"
#define PARAM_MATRUITY_OGIVE_LABEL                "maturity_ogive_label"
#define PARAM_MATURATION                          "maturation"
#define PARAM_MATURATION_RATE                     "maturation_rate"
#define PARAM_MATURE_BIOMASS                      "mature_biomass"
#define PARAM_MAX_AGE                             "max_age"
#define PARAM_MAX_ITER                            "max_iter"
#define PARAM_MAX_THREADS_TO_USE                  "max_threads_to_use"
#define PARAM_MAX_ITERATIONS                      "iterations"
#define PARAM_MAX_LENGTH                          "maximum_length"
#define PARAM_MEAN                                "mean"
#define PARAM_MEAN_LOWER_BOUND                    "mean_lower_bound"
#define PARAM_MEAN_UPPER_BOUND                    "mean_upper_bound"
#define PARAM_MEDIUM                              "medium"
#define PARAM_META                                "meta"
#define PARAM_META_NUMERIC                        "meta_numeric"
#define PARAM_METHOD                              "method"
#define PARAM_METROPOLIS_HASTINGS                 "metropolis_hastings"
#define PARAM_METHOD_OF_REMOVAL                   "method_of_removal"
#define PARAM_MIGRATION_MOVEMENT                  "migration"
#define PARAM_MIN_AGE                             "min_age"
#define PARAM_MINIMUM_LEGAL_LENGTH                "minimum_legel_length"
#define PARAM_MINIMIZER                           "minimiser"
#define PARAM_MODEL                               "model"
#define PARAM_MODEL_ATTRIBUTES                    "model_attributes"
#define PARAM_MORTALITY                           "mortality"
#define PARAM_MORTALITY_CONSTANT_RATE             "mortality_constant_rate"
#define PARAM_MORTALITY_EVENT                     "mortality_event"
#define PARAM_MORTALITY_EVENT_BIOMASS             "mortality_event_biomass"
#define PARAM_MORTALITY_INSTANTANEOUS             "mortality_instantaneous"
#define PARAM_MORTALITY_INITIALISATION_EVENT      "mortality_initialisation_event"
#define PARAM_MORTALITY_INITIALISATION_EVENT_BIOMSS "mortality_initialisation_event_biomass"
#define PARAM_MORTALITY_INSTANTANEOUS_PROCESS     "mortality_instantaneous_process"
#define PARAM_MORTALITY_HOLLING_RATE              "mortality_holling_rate"
#define PARAM_MOVEMENT_BOX_TRANSFER               "movement_box_transfer"
#define PARAM_MOVEMENT_TYPE                       "movement_type"
#define PARAM_MPD                                 "mpd"
#define PARAM_MU                                  "mu"
#define PARAM_MULTINOMIAL                         "multinomial"
#define PARAM_MULTIPLICATIVE                      "multiplicative"
#define PARAM_MULTIPLIER                          "multiplier"
#define PARAM_N                                   "n"
#define PARAM_NAMES                               "names"
#define PARAM_NATAL_HOMING                        "natal_homing"
#define	PARAM_NATURAL_MORTALITY_PROCESS_LABEL	  "natural_mortality_process_label"
#define PARAM_NCOLS                               "ncols"
#define PARAM_NEAREST_NEIGHBOUR                   "nearest_neighbour"
#define PARAM_NONE                                "none"
#define PARAM_NOP                                 "nop"
#define PARAM_NORMAL                              "normal"
#define PARAM_NORMALISED_RESIDUALS                "normalised_residuals"
#define PARAM_NORMAL_BY_STDEV                     "normal_by_stdev"
#define PARAM_NORMAL_LOG                          "normal_log"
#define PARAM_NO_STANDARD_HEADER_REPORT           "no_standard_header_report"
#define PARAM_NUISANCE                            "nuisance"
#define PARAM_NUMBERS                             "numbers"
#define PARAM_NUMBER_OF_GROWTH_EPISODES           "number_of_growth_episodes"
#define PARAM_NUMBER_OF_AGENTS					         "number_of_agents"
#define PARAM_NUMERIC                             "numeric"
#define PARAM_NUMERIC_META                        "numeric_meta"
#define PARAM_NROWS                               "nrows"
#define PARAM_OBJECTIVE                           "objective"
#define PARAM_OBJECTIVE_FUNCTION                  "objective_function"
#define PARAM_OBS                                 "obs"
#define PARAM_OBSERVATION                         "observation"
#define PARAM_OFF_BY_ONE                          "off_by_one"
#define PARAM_OUTPUT_PARAMETERS                   "output_parameters"
#define PARAM_ONE                                 "1 (one)"
#define PARAM_ONE_THOUSAND                        "1000"
#define PARAM_ORIGIN_CELL                         "origin_cell"
#define PARAM_ORTHOGONAL                          "orthogonal"
#define PARAM_OVERWRITE                           "overwrite"
#define PARAM_P1                                  "p1"
#define PARAM_P2                                  "p2"
#define PARAM_PARAMETER                           "parameter"
#define PARAM_PARTITION                           "partition"
#define PARAM_PARTITION_BIOMASS                   "partition_biomass"
#define PARAM_PARTITION_MEAN_WEIGHT               "partition_mean_weight"
#define PARAM_PARTITION_TYPE                      "partition_type"
#define PARAM_PENALTY                             "penalty"
#define PARAM_PEARSONS_RESIDUALS                  "pearsons_residuals"
#define PARAM_PLUS_GROUP                          "plus_group"
#define PARAM_POINT_PERTUBATION_RADIUS            "point_pertubation_radius"
#define PARAM_POPULATION_SIZE                     "population_size"
#define PARAM_PREY_CATEGORIES                     "prey_categories"
#define PARAM_PREY_SUITABILITY_PREDATION          "prey_suitability_predation"
#define PARAM_PREY_SELECTIVITIES_BY_YEAR          "prey_selectivities_by_year"
#define PARAM_PREDATOR_CATEGORIES                 "predator_categories"
#define PARAM_PREY_SELECTIVITIES                  "prey_selectivities"
#define PARAM_PREDATOR_SELECTIVITIES              "predator_selectivities"
#define PARAM_PREDATOR_SELECTIVITIES_BY_YEAR      "predator_selectivities_by_year"
#define PARAM_PREFERENCE_MOVEMENT                 "preference"
#define PARAM_PREFERENCE_FUNCTION                 "preference_function"
#define PARAM_PREFERENCE_FUNCTIONS                "preference_functions"
#define PARAM_PRINT_DEFAULT_REPORTS               "print_default_reports"
#define PARAM_PRINT_LEVEL                         "print_level"
#define PARAM_PRINT_REPORT                        "print_report"
#define PARAM_PRIOR                               "prior"
#define PARAM_PRIOR_APPLIES_TO_TRANSFORM          "prior_applies_to_transform"
#define PARAM_PRIOR_YCS_VALUES                    "prior_standardised_ycs"
#define PARAM_PROCESS                             "process"
#define PARAM_PROCESS_ERROR                       "process_error"
#define PARAM_PROCESS_ERRORS                      "process_errors"
#define PARAM_PROCESS_ERROR_TYPE                  "process_error_type"
#define PARAM_PROCESS_PROPORTION                  "process_proportion"
#define PARAM_PROCESSES                           "processes"
#define PARAM_PROFILE                             "profile"
#define PARAM_PROFILES                            "profiles"
#define PARAM_PROJECT                             "project"
#define PARAM_PROJECTS                            "projects"
#define PARAM_PROCESS_ABUNDANCE                   "process_abundance"
#define PARAM_PROCESS_BIOMASS                     "process_biomass"
#define PARAM_PROCESS_LABEL                       "process_label"
#define PARAM_PROCESS_PROPORTIONS_AT_AGE          "process_proportions_at_age"
#define PARAM_PROCESS_PROPORTIONS_AT_LENGTH       "process_proportions_at_length"
#define PARAM_PROCESS_PROPORTIONS_BY_CATEGORY     "process_proportions_by_category"
#define PARAM_PROCESS_PROPORTIONS_MIGRATING       "process_proportions_migrating"
#define PARAM_PROCESS_REMOVALS_BY_AGE             "process_removals_by_age"
#define PARAM_PROCESS_REMOVALS_BY_LENGTH          "process_removals_by_length"
#define PARAM_PROJECTION_FINAL_YEAR               "projection_final_year"
#define PARAM_PROBABILITY_LAYERS                  "probability_layers"
#define PARAM_PROPORTION                          "proportion"
#define PARAM_PROPORTION_TRHOUGH_MORTALITY        "proportion_through_mortality_block"
#define PARAM_PROPORTIONS                         "proportions"
#define PARAM_PROPORTIONS_AT_AGE                  "proportions_at_age"
#define PARAM_PROPORTIONS_AT_LENGTH               "proportions_at_length"
#define PARAM_PROPORTIONS_MATURE_BY_AGE           "proportions_mature_by_age"
#define PARAM_PROPORTION_MALE                     "proportion_male"
#define PARAM_PROPORTIONS_BY_CATEGORY             "proportions_by_category"
#define PARAM_PROPORTION_TIME_STEP                "proportion_time-step"
#define PARAM_PROPOSAL_DISTRIBUTION               "proposal_distribution"
#define PARAM_PSEUDO                              "none"
#define PARAM_PSI                                 "psi"
#define PARAM_Q                                   "q"
#define PARAM_R                                   "r"
#define PARAM_RECAPTURED                          "recaptured"
#define PARAM_R0                                  "r0"
#define PARAM_R0_LAYER                            "r0_layer"
#define PARAM_RANDOM_NUMBER_SEED                  "random_number_seed"
#define PARAM_RANDOMDRAW                          "random_draw"
#define PARAM_RANDOMWALK                          "random_walk"
#define PARAM_RATE                                "rate"
#define PARAM_RATES                               "rates"
#define PARAM_RATIO                               "ratio"
#define PARAM_RECRUITMENT                         "recruitment"
#define PARAM_RECRUITMENT_TIME                    "recruitment_before_ageing"
#define PARAM_RECRUITMENT_BEVERTON_HOLT           "recruitment_beverton_holt"
#define PARAM_RECRUITMENT_BEVERTON_HOLT_WITH_DEVIATIONS "recruitment_beverton_holt_with_deviations"
#define PARAM_RECRUITMENT_CONSTANT                "recruitment_constant"
#define PARAM_RECRUITMENT_LABEL                   "recruitment_label"
#define PARAM_RECRUITMENT_VALUES                  "recruitment_values"
#define PARAM_RECRUITMENT_LAYER_LABEL             "recruitment_layer_label"
#define PARAM_REPORT                              "report"
#define PARAM_RESCALE                             "rescale"
#define PARAM_RETAPE                              "retape"
#define PARAM_RHO                                 "rho"
#define PARAM_ROW                                 "row"
#define PARAM_ROW_INDEX                           "row_index"
#define PARAM_RUN_LENGTH                          "run_length"
#define PARAM_RUN_MODE                            "run_mode"
#define PARAM_S                                   "s"
#define PARAM_SAME                                "same"
#define PARAM_SAMPLE                              "sample"
#define PARAM_SB                                  "sb"
#define PARAM_SCALING_YEARS                       "scaling_years"
#define PARAM_SCANNED                             "scanned"
#define PARAM_SCORES_INDEX                        "scores_index"
#define PARAM_SCHNUTE                             "schnute"
#define PARAM_SECTION                             "command"
#define PARAM_SECOND_ESTIMATE                     "second_estimate"
#define PARAM_SECOND_PARAMETER                    "second_parameter"
#define PARAM_SELECTIVITIES                       "selectivities"
#define PARAM_SELECTIVITY                         "selectivity"
#define PARAM_SELECTIVITY_LABEL                   "selectivity_label"
#define PARAM_SEQUENTIALLY_ADD                    "sequentially_add"
#define PARAM_SEXED                               "sexed"
#define PARAM_SEPERATE_BY_SEX                     "seperate_by_sex"
#define PARAM_SHIFT_A                             "shift_parameter"
#define PARAM_SIGMA                               "sigma"
#define PARAM_SIGMA_L                             "sigma_l"
#define PARAM_SIGMA_R                             "sigma_r"
#define PARAM_SIGMA_MIN                           "sigma_min"
#define PARAM_SIMPLEX                             "simplex"
#define PARAM_SIMULATED_OBSERVATION               "simulated_observation"
#define PARAM_SIMULATION_LIKELIHOOD               "simulation_likelihood"
#define PARAM_SINGLE_STEP                         "single_step"
#define PARAM_SINGLE                              "single"
#define PARAM_SINK_LAYER                          "sink_layer"
#define PARAM_SIZES                               "sizes"
#define PARAM_SLOPE                               "slope"
#define PARAM_SKIP_CONFIG_FILE                    "skip_config_file"
#define PARAM_SOLVER                              "solver"
#define PARAM_SOURCE_LAYER                        "source_layer"
#define PARAM_SPATIAL_MAP                         "spatial_map"
#define PARAM_SSB_LAYER                           "ssb_layer"
#define PARAM_SSB_VALUES                          "ssb_values"
#define PARAM_STANDARDISE_YCS_YEARS               "standardise_ycs_years"
#define PARAM_STANDARDISE_YCS                     "standardise_ycs"
#define PARAM_START                               "start"
#define PARAM_START_YEAR                          "start_year"
#define PARAM_STATE                               "state"
#define PARAM_SQUARE_ROOT                         "sqrt"
#define PARAM_SQUARE_UNIFORM                      "square_uniform"
#define PARAM_SSB                                 "ssb"
#define PARAM_SSB_OFFSET                          "ssb_offset"
#define PARAM_SSB_LAYER_LABEL                     "ssb_layer_label"
#define PARAM_STATE_CATEGORY_BY_AGE               "state_category_by_age"
#define PARAM_STEP_SIZE                           "step_size"
#define PARAM_STEPS                               "steps"
#define PARAM_STEEPNESS                           "steepness"
#define PARAM_STRING                              "categorical"
#define PARAM_SUM_TO_ONE                          "sum_to_one"
#define PARAM_SUMMARISE_AGENTS                    "summarise_agents"
#define PARAM_T                                   "t"
#define PARAM_T0                                  "t0"
#define PARAM_TABLE                               "table"
#define PARAM_TAG                                 "tag"
#define PARAM_TAGGED_CATEGORIES                   "tagged_categories"
#define PARAM_TAGGED_SELECTIVITIES                "tagged_selectivities"
#define PARAM_TAG_BY_AGE                          "tag_by_age"
#define PARAM_TAG_BY_LENGTH                       "tag_by_length"
#define PARAM_TAG_LOSS_RATE                       "tag_loss_rate"
#define PARAM_TAG_LOSS_YEARS                      "tag_loss_years"
#define PARAM_TAG_LOSS_TYPE                       "tag_loss_type"
#define PARAM_TAG_LOSS                            "tag_loss"
#define PARAM_TAG_RECAPTURE_BY_AGE                "tag_recapture_by_age"
#define PARAM_TAG_RECAPTURE_BY_LENGTH             "tag_recapture_by_length"
#define PARAM_TAU1                                "tau1"
#define PARAM_TAU2                                "tau2"
#define PARAM_TARGET_CATEGORIES                   "categories2"
#define PARAM_TARGET_SELECTIVITIES                "selectivities2"
#define PARAM_TERMINAL_YEAR                       "terminal_year"
#define PARAM_THRESHOLD                           "threshold"
#define PARAM_THRESHOLD_BIOMASS                   "threshold_biomass"
#define PARAM_THETA_ONE                           "theta1"
#define PARAM_THETA_TWO                           "theta2"
#define PARAM_TIME_STEP                           "time_step"
#define PARAM_TIME_STEP_MEASUREMENTS_WERE_MADE    "time_step_measurements_were_made"
#define PARAM_TIME_STEP_PROPORTION                "time_step_proportion"
#define PARAM_TIME_STEP_PROPORTIONS               "time_step_proportions"
#define PARAM_TIME_STEP_PROPORTION_METHOD         "time_step_proportion_method"
#define PARAM_TIME_STEP_RATIO                     "time_step_ratio"
#define PARAM_TIME_STEPS                          "time_steps"
#define PARAM_TIME_VARYING                        "time_varying"
#define PARAM_TO                                  "to"
#define PARAM_TOL                                 "tol"
#define PARAM_TOLERANCE                           "tolerance"
#define PARAM_TONNES                              "tonnes"
#define PARAM_TOTAL_CATEGORIES                    "total_categories"
#define PARAM_TOTAL_SCORE                         "total_score"
#define PARAM_TOTAL                               "total"
#define PARAM_TRACE                               "trace"
#define PARAM_TRANSFORM_WITH_JACOBIAN             "transform_with_jacobian"
#define PARAM_TRANSFORMATION                      "transformation"
#define PARAM_TRANSFORMATION_TYPE                 "transformation_type"
#define PARAM_TRANSITION                          "transition"
#define PARAM_TRANSITION_CATEGORY                 "transition_category"
#define PARAM_TRANSITION_CATEGORY_BY_AGE          "transition_category_by_age"
#define PARAM_TRUE_YCS_VALUES                     "true_ycs_values"
#define PARAM_TYPE                                "type"
#define PARAM_U_MAX                               "u_max"
#define PARAM_UOBS_F                              "Fishing_pressure"
#define PARAM_UNITS                               "units"
#define PARAM_UNIFORM                             "uniform"
#define PARAM_UNIFORM_LOG                         "uniform_log"
#define PARAM_UPPER_BOUND                         "upper_bound"
#define PARAM_UPDATE_GROWTH_PARAMETERS            "update_growth_parameters"
#define PARAM_UPDATE_MORTALITY_PARAMETERS         "update_mortality_parameters"
#define PARAM_USER_DEFINED                        "user_defined"
#define PARAM_V                                   "v"
#define PARAM_VECTOR_AVERAGE                      "vector_average"
#define PARAM_VECTOR_SMOOTHING                    "vector_smoothing"
#define PARAM_VON_BERTALANFFY                     "von_bertalanffy"
#define PARAM_VALUE                               "value"
#define PARAM_VALUES                              "values"
#define PARAM_WIDTH                               "width"
#define PARAM_WEIGHTS                             "weights"
#define PARAM_WEIGHTED_PRODUCT                    "weighted_product"
#define PARAM_WEIGHTED_SUM                        "weighted_sum"
#define PARAM_WRITE_MODE                          "write_mode"
#define PARAM_WORLD_AGE_FREQUENCY                 "world_age_frequency"
#define PARAM_X0                                  "x0"
#define PARAM_X                                   "x"
#define PARAM_X1                                  "x1"
#define PARAM_X2                                  "x2"
#define PARAM_Y0                                  "y0"
#define PARAM_Y1                                  "y1"
#define PARAM_Y2                                  "y2"
#define PARAM_YCS_VALUES                          "ycs_values"
#define PARAM_YCS_YEARS                           "ycs_years"
#define PARAM_YEAR                                "year"
#define PARAM_YEARS                               "years"
#define PARAM_ZERO                                "0 (zero)"

#endif /* ENGLISH_H_ */
