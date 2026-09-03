#pragma once

/*
 * BDTReweighter_json.hh
 *
 * Header-only C++17 evaluator for hep_ml GBReweighter models exported by
 * json_reweighter.py.
 *
 * JSON dependency:
 *     Boost.JSON (tested against the Boost 1.82 API)
 *
 * Link the final program with Boost.JSON:
 *
 *     g++ -std=c++17 -O3 -march=native -ffp-contract=off program.cc \
 *         -I"$BOOST_INC" -L"$BOOST_LIB" \
 *         -Wl,-rpath,"$BOOST_LIB" -lboost_json -o program
 *
 * This header intentionally does NOT include <boost/json/src.hpp>.
 *
 * The input feature vector is positional. Its entries must have the same
 * order as the columns used to train the Python reweighter.
 *
 * The JSON is parsed only in the constructor. Trees are then stored in one
 * compact typed node pool with global child indices, so single-event
 * prediction performs no JSON access, allocation, locking, or event-level
 * parallelism. Supplying float32 features selects the lowest-overhead path;
 * double features remain supported and are converted at each visited split.
 */

#include <boost/json.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "systematicstools/utility/string_parsers.hh"

namespace BDTReweight {

namespace detail {

inline constexpr const char* kModelFormat =
    "hep_ml_gb_reweighter";

inline constexpr std::int64_t kFormatVersion = 1;

[[nodiscard]] inline std::runtime_error ModelError(
    const std::string& message
) {
    return std::runtime_error(
        "JSONReweighter: " + message
    );
}

[[nodiscard]] inline std::string TreeContext(
    std::size_t tree_index,
    std::size_t node_index
) {
    std::ostringstream stream;
    stream
        << "tree " << tree_index
        << ", node " << node_index;
    return stream.str();
}

[[nodiscard]] inline const boost::json::value&
GetRequiredValue(
    const boost::json::object& object,
    boost::json::string_view key,
    const std::string& context
) {
    const boost::json::value* value =
        object.if_contains(key);

    if (value == nullptr) {
        throw ModelError(
            context + " is missing field '"
            + std::string(key.data(), key.size()) + "'."
        );
    }

    return *value;
}

[[nodiscard]] inline const boost::json::object&
AsObject(
    const boost::json::value& value,
    const std::string& context
) {
    if (!value.is_object()) {
        throw ModelError(
            context + " must be a JSON object."
        );
    }

    return value.as_object();
}

[[nodiscard]] inline const boost::json::array&
AsArray(
    const boost::json::value& value,
    const std::string& context
) {
    if (!value.is_array()) {
        throw ModelError(
            context + " must be a JSON array."
        );
    }

    return value.as_array();
}

[[nodiscard]] inline std::string AsString(
    const boost::json::value& value,
    const std::string& context
) {
    if (!value.is_string()) {
        throw ModelError(
            context + " must be a JSON string."
        );
    }

    const boost::json::string& text =
        value.as_string();

    return std::string(text.data(), text.size());
}

[[nodiscard]] inline bool AsBool(
    const boost::json::value& value,
    const std::string& context
) {
    if (!value.is_bool()) {
        throw ModelError(
            context + " must be a JSON boolean."
        );
    }

    return value.as_bool();
}

[[nodiscard]] inline std::int64_t AsInt64(
    const boost::json::value& value,
    const std::string& context
) {
    if (value.is_int64()) {
        return value.as_int64();
    }

    if (value.is_uint64()) {
        const std::uint64_t number =
            value.as_uint64();

        if (
            number
            <= static_cast<std::uint64_t>(
                std::numeric_limits<std::int64_t>::max()
            )
        ) {
            return static_cast<std::int64_t>(number);
        }
    }

    throw ModelError(
        context + " must be an integer representable "
        "as int64."
    );
}

[[nodiscard]] inline double AsDouble(
    const boost::json::value& value,
    const std::string& context
) {
    double result = 0.0;

    if (value.is_double()) {
        result = value.as_double();
    } else if (value.is_int64()) {
        result = static_cast<double>(
            value.as_int64()
        );
    } else if (value.is_uint64()) {
        result = static_cast<double>(
            value.as_uint64()
        );
    } else {
        throw ModelError(
            context + " must be a JSON number."
        );
    }

    if (!std::isfinite(result)) {
        throw ModelError(
            context + " must be finite."
        );
    }

    return result;
}

[[nodiscard]] inline std::string GetString(
    const boost::json::object& object,
    boost::json::string_view key,
    const std::string& context
) {
    return AsString(
        GetRequiredValue(object, key, context),
        context + " field '"
            + std::string(key.data(), key.size()) + "'"
    );
}

[[nodiscard]] inline bool GetBool(
    const boost::json::object& object,
    boost::json::string_view key,
    const std::string& context
) {
    return AsBool(
        GetRequiredValue(object, key, context),
        context + " field '"
            + std::string(key.data(), key.size()) + "'"
    );
}

[[nodiscard]] inline std::int64_t GetInt64(
    const boost::json::object& object,
    boost::json::string_view key,
    const std::string& context
) {
    return AsInt64(
        GetRequiredValue(object, key, context),
        context + " field '"
            + std::string(key.data(), key.size()) + "'"
    );
}

[[nodiscard]] inline double GetDouble(
    const boost::json::object& object,
    boost::json::string_view key,
    const std::string& context
) {
    return AsDouble(
        GetRequiredValue(object, key, context),
        context + " field '"
            + std::string(key.data(), key.size()) + "'"
    );
}

[[nodiscard]] inline boost::json::value ParseFile(
    const std::filesystem::path& filepath
) {
    std::ifstream input( systtools::expand_env_vars(filepath) );

    if (!input) {
        throw ModelError(
            "could not open model file '"
            + filepath.string() + "'."
        );
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();

    if (input.bad()) {
        throw ModelError(
            "failed while reading model file '"
            + filepath.string() + "'."
        );
    }

    try {
        return boost::json::parse(buffer.str());
    } catch (const std::exception& error) {
        throw ModelError(
            "failed to parse '" + filepath.string()
            + "': " + error.what()
        );
    }
}

}  // namespace detail


/**
 * @brief Standalone evaluator for a hep_ml GBReweighter exported to JSON.
 *
 * For every split, the selected feature value is converted to IEEE-754
 * float32 while the stored threshold remains IEEE-754 float64. This is the
 * traversal behavior verified against scikit-learn for the exported model.
 * The compiled model is immutable after construction, so concurrent calls on
 * distinct events are safe if the surrounding application later adds its own
 * event-level parallelism.
 */
class JSONReweighter {
public:
    /**
     * @brief Load and validate an exported JSON model.
     *
     * @param filepath Path to a model produced by
     *        export_reweighter_json() in json_reweighter.py.
     *
     * @throws std::runtime_error for file, parsing, or model-schema errors.
     */
    explicit JSONReweighter(
        const std::filesystem::path& filepath
    );

    [[nodiscard]] std::size_t NFeatures() const noexcept;
    [[nodiscard]] std::size_t NTrees() const noexcept;
    [[nodiscard]] double LearningRate() const noexcept;
    [[nodiscard]] double InitialStep() const noexcept;

    /**
     * @brief Evaluate the ensemble score before exponentiation.
     */
    [[nodiscard]] double DecisionFunction(
        const std::vector<double>& features
    ) const;

    /**
     * @brief Evaluate the ensemble score from a raw feature buffer.
     */
    [[nodiscard]] double DecisionFunction(
        const double* features,
        std::size_t feature_count
    ) const;

    /**
     * @brief Evaluate from already-converted float32 features.
     *
     * This is the fastest scalar path because no conversion is needed while
     * traversing the trees.
     */
    [[nodiscard]] double DecisionFunction(
        const std::vector<float>& features
    ) const;

    [[nodiscard]] double DecisionFunction(
        const float* features,
        std::size_t feature_count
    ) const;

    /**
     * @brief Predict one event weight.
     *
     * With the default original_weight=1, this returns the learned
     * multiplicative reweight.
     */
    [[nodiscard]] double PredictWeight(
        const std::vector<double>& features,
        double original_weight = 1.0
    ) const;

    /**
     * @brief Predict one event weight from a raw feature buffer.
     */
    [[nodiscard]] double PredictWeight(
        const double* features,
        std::size_t feature_count,
        double original_weight = 1.0
    ) const;

    /**
     * @brief Predict from already-converted float32 features.
     */
    [[nodiscard]] double PredictWeight(
        const std::vector<float>& features,
        double original_weight = 1.0
    ) const;

    [[nodiscard]] double PredictWeight(
        const float* features,
        std::size_t feature_count,
        double original_weight = 1.0
    ) const;

    /**
     * @brief Predict multipliers for a flattened row-major event matrix.
     *
     * events must contain n_events * NFeatures() values.
     */
    [[nodiscard]] std::vector<double> PredictWeights(
        const std::vector<double>& events,
        std::size_t n_events
    ) const;

    /**
     * @brief Predict updated event weights for a flattened row-major matrix.
     */
    [[nodiscard]] std::vector<double> PredictWeights(
        const std::vector<double>& events,
        std::size_t n_events,
        const std::vector<double>& original_weights
    ) const;

private:
    using NodeIndex = std::uint32_t;

    static constexpr std::int32_t kLeafFeature = -1;

    // Compact inference representation. For an internal node,
    // split_or_leaf_value is its float64 threshold. For a leaf, it is the
    // leaf contribution and feature == kLeafFeature.
    struct Node {
        double split_or_leaf_value = 0.0;
        NodeIndex left = 0;
        NodeIndex right = 0;
        std::int32_t feature = kLeafFeature;
    };

    struct Tree {
        std::vector<Node> nodes;
    };

    template <typename Feature>
    [[nodiscard]] double EvaluateTree(
        NodeIndex root,
        const Feature* features
    ) const noexcept;

    template <typename Feature>
    [[nodiscard]] double DecisionFunctionUnchecked(
        const Feature* features
    ) const noexcept;

    static void ValidateTreeGraph(
        const Tree& tree,
        std::size_t tree_index
    );

    template <typename Feature>
    void ValidateFeatures(
        const Feature* features,
        std::size_t feature_count
    ) const;

    void ValidateFlatMatrix(
        const std::vector<double>& events,
        std::size_t n_events
    ) const;

    std::size_t n_features_ = 0;
    double learning_rate_ = 0.0;
    double initial_step_ = 0.0;
    std::vector<Node> nodes_;
    std::vector<NodeIndex> tree_roots_;
};


inline JSONReweighter::JSONReweighter(
    const std::filesystem::path& filepath
) {
    static_assert(
        sizeof(float) == 4
            && std::numeric_limits<float>::is_iec559,
        "JSONReweighter requires IEEE-754 32-bit float."
    );

    static_assert(
        sizeof(double) == 8
            && std::numeric_limits<double>::is_iec559,
        "JSONReweighter requires IEEE-754 64-bit double."
    );

    const boost::json::value document =
        detail::ParseFile(filepath);

    const std::string top_level_context =
        "top-level model";

    const boost::json::object& model =
        detail::AsObject(document, top_level_context);

    const std::string format =
        detail::GetString(model, "format", "model");

    if (format != detail::kModelFormat) {
        throw detail::ModelError(
            "unsupported model format '" + format + "'."
        );
    }

    const std::int64_t format_version =
        detail::GetInt64(
            model,
            "format_version",
            "model"
        );

    if (format_version != detail::kFormatVersion) {
        throw detail::ModelError(
            "unsupported format version "
            + std::to_string(format_version)
            + "; expected "
            + std::to_string(detail::kFormatVersion)
            + "."
        );
    }

    const std::int64_t n_features =
        detail::GetInt64(
            model,
            "n_features",
            "model"
        );

    if (n_features <= 0) {
        throw detail::ModelError(
            "n_features must be positive."
        );
    }

    if (
        static_cast<std::uint64_t>(n_features)
        > static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max()
        )
    ) {
        throw detail::ModelError(
            "n_features is too large for this platform."
        );
    }

    if (
        n_features
        > static_cast<std::int64_t>(
            std::numeric_limits<int>::max()
        )
    ) {
        throw detail::ModelError(
            "n_features is too large for this evaluator."
        );
    }

    n_features_ =
        static_cast<std::size_t>(n_features);

    const std::int64_t n_trees =
        detail::GetInt64(
            model,
            "n_trees",
            "model"
        );

    if (n_trees < 0) {
        throw detail::ModelError(
            "n_trees must be nonnegative."
        );
    }

    learning_rate_ =
        detail::GetDouble(
            model,
            "learning_rate",
            "model"
        );

    initial_step_ =
        detail::GetDouble(
            model,
            "initial_step",
            "model"
        );

    if (
        const boost::json::value* input_precision =
            model.if_contains("input_precision")
    ) {
        const std::string precision =
            detail::AsString(
                *input_precision,
                "model field 'input_precision'"
            );

        if (precision != "float32") {
            throw detail::ModelError(
                "unsupported input_precision '"
                + precision + "'."
            );
        }
    }

    if (
        const boost::json::value* threshold_precision =
            model.if_contains("threshold_precision")
    ) {
        const std::string precision =
            detail::AsString(
                *threshold_precision,
                "model field 'threshold_precision'"
            );

        if (precision != "float64") {
            throw detail::ModelError(
                "unsupported threshold_precision '"
                + precision + "'."
            );
        }
    }

    if (
        const boost::json::value* non_finite_policy =
            model.if_contains("non_finite_policy")
    ) {
        const std::string policy =
            detail::AsString(
                *non_finite_policy,
                "model field 'non_finite_policy'"
            );

        if (policy != "reject") {
            throw detail::ModelError(
                "unsupported non_finite_policy '"
                + policy + "'."
            );
        }
    }

    const std::string trees_context =
        "model field 'trees'";
    const std::string model_context = "model";

    const boost::json::value& trees_value =
        detail::GetRequiredValue(
            model,
            "trees",
            model_context
        );

    const boost::json::array& json_trees =
        detail::AsArray(
            trees_value,
            trees_context
        );

    if (
        static_cast<std::uint64_t>(n_trees)
        != static_cast<std::uint64_t>(
            json_trees.size()
        )
    ) {
        throw detail::ModelError(
            "n_trees is " + std::to_string(n_trees)
            + ", but the JSON contains "
            + std::to_string(json_trees.size())
            + " trees."
        );
    }

    tree_roots_.reserve(json_trees.size());

    for (
        std::size_t tree_index = 0;
        tree_index < json_trees.size();
        ++tree_index
    ) {
        const std::string tree_context =
            "tree " + std::to_string(tree_index);

        const boost::json::object& json_tree =
            detail::AsObject(
                json_trees.at(tree_index),
                tree_context
            );

        const std::int64_t node_count =
            detail::GetInt64(
                json_tree,
                "node_count",
                tree_context
            );

        if (node_count <= 0) {
            throw detail::ModelError(
                tree_context
                + " has nonpositive node_count."
            );
        }

        if (
            node_count
            > static_cast<std::int64_t>(
                std::numeric_limits<int>::max()
            )
        ) {
            throw detail::ModelError(
                tree_context
                + " has too many nodes for this evaluator."
            );
        }

        const std::string nodes_context =
            tree_context + " field 'nodes'";

        const boost::json::array& json_nodes =
            detail::AsArray(
                detail::GetRequiredValue(
                    json_tree,
                    "nodes",
                    tree_context
                ),
                nodes_context
            );

        if (
            static_cast<std::uint64_t>(node_count)
            != static_cast<std::uint64_t>(
                json_nodes.size()
            )
        ) {
            throw detail::ModelError(
                tree_context + " declares "
                + std::to_string(node_count)
                + " nodes, but contains "
                + std::to_string(json_nodes.size())
                + "."
            );
        }

        Tree tree;
        tree.nodes.reserve(json_nodes.size());

        for (
            std::size_t node_index = 0;
            node_index < json_nodes.size();
            ++node_index
        ) {
            const std::string node_context =
                detail::TreeContext(
                    tree_index,
                    node_index
                );

            const boost::json::object& json_node =
                detail::AsObject(
                    json_nodes.at(node_index),
                    node_context
                );

            Node node;
            const bool is_leaf = detail::GetBool(
                json_node,
                "is_leaf",
                node_context
            );

            if (is_leaf) {
                node.split_or_leaf_value =
                    detail::GetDouble(
                        json_node,
                        "value",
                        node_context
                    );
            } else {
                const std::int64_t feature =
                    detail::GetInt64(
                        json_node,
                        "feature",
                        node_context
                    );

                const std::int64_t left =
                    detail::GetInt64(
                        json_node,
                        "left",
                        node_context
                    );

                const std::int64_t right =
                    detail::GetInt64(
                        json_node,
                        "right",
                        node_context
                    );

                node.split_or_leaf_value =
                    detail::GetDouble(
                        json_node,
                        "threshold",
                        node_context
                    );

                if (
                    feature < 0
                    || static_cast<std::uint64_t>(feature)
                        >= static_cast<std::uint64_t>(
                            n_features_
                        )
                ) {
                    throw detail::ModelError(
                        node_context
                        + " has feature index "
                        + std::to_string(feature)
                        + " outside [0, "
                        + std::to_string(n_features_)
                        + ")."
                    );
                }

                if (
                    left < 0
                    || static_cast<std::uint64_t>(left)
                        >= static_cast<std::uint64_t>(
                            json_nodes.size()
                        )
                ) {
                    throw detail::ModelError(
                        node_context
                        + " has invalid left child "
                        + std::to_string(left) + "."
                    );
                }

                if (
                    right < 0
                    || static_cast<std::uint64_t>(right)
                        >= static_cast<std::uint64_t>(
                            json_nodes.size()
                        )
                ) {
                    throw detail::ModelError(
                        node_context
                        + " has invalid right child "
                        + std::to_string(right) + "."
                    );
                }

                node.feature =
                    static_cast<std::int32_t>(feature);
                node.left =
                    static_cast<NodeIndex>(left);
                node.right =
                    static_cast<NodeIndex>(right);
            }

            tree.nodes.push_back(node);
        }

        ValidateTreeGraph(tree, tree_index);

        const std::size_t max_node_count =
            static_cast<std::size_t>(
                std::numeric_limits<NodeIndex>::max()
            );

        if (
            tree.nodes.size()
            > max_node_count - nodes_.size()
        ) {
            throw detail::ModelError(
                "the combined model has too many nodes "
                "for this evaluator."
            );
        }

        const NodeIndex node_offset =
            static_cast<NodeIndex>(nodes_.size());

        tree_roots_.push_back(node_offset);

        for (Node node : tree.nodes) {
            if (node.feature != kLeafFeature) {
                node.left = static_cast<NodeIndex>(
                    node.left + node_offset
                );
                node.right = static_cast<NodeIndex>(
                    node.right + node_offset
                );
            }

            nodes_.push_back(node);
        }
    }
}


inline std::size_t
JSONReweighter::NFeatures() const noexcept {
    return n_features_;
}


inline std::size_t
JSONReweighter::NTrees() const noexcept {
    return tree_roots_.size();
}


inline double
JSONReweighter::LearningRate() const noexcept {
    return learning_rate_;
}


inline double
JSONReweighter::InitialStep() const noexcept {
    return initial_step_;
}


inline void JSONReweighter::ValidateTreeGraph(
    const Tree& tree,
    std::size_t tree_index
) {
    const std::size_t node_count =
        tree.nodes.size();

    std::vector<int> parent_count(
        node_count,
        0
    );

    for (
        std::size_t node_index = 0;
        node_index < node_count;
        ++node_index
    ) {
        const Node& node =
            tree.nodes[node_index];

        if (node.feature == kLeafFeature) {
            continue;
        }

        if (node.left == node.right) {
            throw detail::ModelError(
                detail::TreeContext(
                    tree_index,
                    node_index
                )
                + " has identical left and right children."
            );
        }

        ++parent_count.at(
            static_cast<std::size_t>(node.left)
        );

        ++parent_count.at(
            static_cast<std::size_t>(node.right)
        );
    }

    if (parent_count.at(0) != 0) {
        throw detail::ModelError(
            "tree " + std::to_string(tree_index)
            + " has a root node with a parent."
        );
    }

    for (
        std::size_t node_index = 1;
        node_index < node_count;
        ++node_index
    ) {
        if (parent_count[node_index] != 1) {
            throw detail::ModelError(
                detail::TreeContext(
                    tree_index,
                    node_index
                )
                + " must have exactly one parent."
            );
        }
    }

    std::vector<int> state(node_count, 0);

    std::function<void(std::size_t)> visit =
        [&](std::size_t node_index) {
            if (state[node_index] == 1) {
                throw detail::ModelError(
                    "tree "
                    + std::to_string(tree_index)
                    + " contains a cycle."
                );
            }

            if (state[node_index] == 2) {
                return;
            }

            state[node_index] = 1;

            const Node& node =
                tree.nodes[node_index];

            if (node.feature != kLeafFeature) {
                visit(
                    static_cast<std::size_t>(
                        node.left
                    )
                );

                visit(
                    static_cast<std::size_t>(
                        node.right
                    )
                );
            }

            state[node_index] = 2;
        };

    visit(0);

    for (
        std::size_t node_index = 0;
        node_index < node_count;
        ++node_index
    ) {
        if (state[node_index] != 2) {
            throw detail::ModelError(
                detail::TreeContext(
                    tree_index,
                    node_index
                )
                + " is unreachable from the root."
            );
        }
    }
}


template <typename Feature>
inline double JSONReweighter::EvaluateTree(
    NodeIndex root,
    const Feature* features
) const noexcept {
    NodeIndex node_index = root;

    while (true) {
        const Node& node =
            nodes_[node_index];

        if (node.feature == kLeafFeature) {
            return node.split_or_leaf_value;
        }

        const double raw_value =
            features[
                static_cast<std::size_t>(
                    node.feature
                )
            ];

        // Precision rule verified against sklearn:
        // feature -> float32, threshold stays float64.
        const float feature_value =
            static_cast<float>(raw_value);

        node_index =
            static_cast<double>(feature_value)
                <= node.split_or_leaf_value
            ? node.left
            : node.right;
    }
}


template <typename Feature>
inline double JSONReweighter::DecisionFunctionUnchecked(
    const Feature* features
) const noexcept {
    double score = initial_step_;

    for (const NodeIndex root : tree_roots_) {
        const double contribution =
            EvaluateTree(root, features);

        const double increment =
            learning_rate_ * contribution;

        // Keep the operation order explicit to match the Python evaluator.
        score = score + increment;
    }

    return score;
}


template <typename Feature>
inline void JSONReweighter::ValidateFeatures(
    const Feature* features,
    std::size_t feature_count
) const {
    if (
        features == nullptr
        && feature_count != 0
    ) {
        throw std::invalid_argument(
            "JSONReweighter: feature pointer is null."
        );
    }

    if (feature_count != n_features_) {
        throw std::invalid_argument(
            "JSONReweighter: expected "
            + std::to_string(n_features_)
            + " features, received "
            + std::to_string(feature_count)
            + "."
        );
    }

    for (
        std::size_t feature_index = 0;
        feature_index < feature_count;
        ++feature_index
    ) {
        if (!std::isfinite(features[feature_index])) {
            throw std::invalid_argument(
                "JSONReweighter: non-finite input "
                "at feature "
                + std::to_string(feature_index)
                + "."
            );
        }
    }
}


inline void JSONReweighter::ValidateFlatMatrix(
    const std::vector<double>& events,
    std::size_t n_events
) const {
    if (
        n_events
        > std::numeric_limits<std::size_t>::max()
            / n_features_
    ) {
        throw std::invalid_argument(
            "JSONReweighter: event-matrix "
            "dimensions overflow."
        );
    }

    const std::size_t expected_size =
        n_events * n_features_;

    if (events.size() != expected_size) {
        throw std::invalid_argument(
            "JSONReweighter: expected a flattened "
            "matrix with "
            + std::to_string(expected_size)
            + " values, received "
            + std::to_string(events.size())
            + "."
        );
    }

    for (
        std::size_t value_index = 0;
        value_index < events.size();
        ++value_index
    ) {
        if (!std::isfinite(events[value_index])) {
            const std::size_t event_index =
                value_index / n_features_;

            const std::size_t feature_index =
                value_index % n_features_;

            throw std::invalid_argument(
                "JSONReweighter: non-finite input "
                "at event "
                + std::to_string(event_index)
                + ", feature "
                + std::to_string(feature_index)
                + "."
            );
        }
    }
}


inline double JSONReweighter::DecisionFunction(
    const std::vector<double>& features
) const {
    return DecisionFunction(
        features.data(),
        features.size()
    );
}


inline double JSONReweighter::DecisionFunction(
    const double* features,
    std::size_t feature_count
) const {
    ValidateFeatures(features, feature_count);
    return DecisionFunctionUnchecked(features);
}


inline double JSONReweighter::DecisionFunction(
    const std::vector<float>& features
) const {
    return DecisionFunction(
        features.data(),
        features.size()
    );
}


inline double JSONReweighter::DecisionFunction(
    const float* features,
    std::size_t feature_count
) const {
    ValidateFeatures(features, feature_count);
    return DecisionFunctionUnchecked(features);
}


inline double JSONReweighter::PredictWeight(
    const std::vector<double>& features,
    double original_weight
) const {
    return PredictWeight(
        features.data(),
        features.size(),
        original_weight
    );
}


inline double JSONReweighter::PredictWeight(
    const double* features,
    std::size_t feature_count,
    double original_weight
) const {
    if (!std::isfinite(original_weight)) {
        throw std::invalid_argument(
            "JSONReweighter: original_weight "
            "must be finite."
        );
    }

    ValidateFeatures(features, feature_count);

    return original_weight
        * std::exp(
            DecisionFunctionUnchecked(features)
        );
}


inline double JSONReweighter::PredictWeight(
    const std::vector<float>& features,
    double original_weight
) const {
    return PredictWeight(
        features.data(),
        features.size(),
        original_weight
    );
}


inline double JSONReweighter::PredictWeight(
    const float* features,
    std::size_t feature_count,
    double original_weight
) const {
    if (!std::isfinite(original_weight)) {
        throw std::invalid_argument(
            "JSONReweighter: original_weight "
            "must be finite."
        );
    }

    ValidateFeatures(features, feature_count);

    return original_weight
        * std::exp(
            DecisionFunctionUnchecked(features)
        );
}


inline std::vector<double>
JSONReweighter::PredictWeights(
    const std::vector<double>& events,
    std::size_t n_events
) const {
    ValidateFlatMatrix(events, n_events);

    std::vector<double> predictions(n_events);

    for (
        std::size_t event_index = 0;
        event_index < n_events;
        ++event_index
    ) {
        const double* event =
            events.data()
            + event_index * n_features_;

        predictions[event_index] = std::exp(
            DecisionFunctionUnchecked(event)
        );
    }

    return predictions;
}


inline std::vector<double>
JSONReweighter::PredictWeights(
    const std::vector<double>& events,
    std::size_t n_events,
    const std::vector<double>& original_weights
) const {
    ValidateFlatMatrix(events, n_events);

    if (original_weights.size() != n_events) {
        throw std::invalid_argument(
            "JSONReweighter: expected "
            + std::to_string(n_events)
            + " original weights, received "
            + std::to_string(
                original_weights.size()
            )
            + "."
        );
    }

    for (
        std::size_t event_index = 0;
        event_index < n_events;
        ++event_index
    ) {
        const double original_weight =
            original_weights[event_index];

        if (!std::isfinite(original_weight)) {
            throw std::invalid_argument(
                "JSONReweighter: non-finite "
                "original weight at event "
                + std::to_string(event_index)
                + "."
            );
        }

    }

    std::vector<double> predictions(n_events);

    for (
        std::size_t event_index = 0;
        event_index < n_events;
        ++event_index
    ) {
        const double* event =
            events.data()
            + event_index * n_features_;

        predictions[event_index] =
            original_weights[event_index]
            * std::exp(
                DecisionFunctionUnchecked(event)
            );
    }

    return predictions;
}

}  // namespace BDTReweight

