#pragma once

#include <mitsuba/core/fwd.h>
#include <algorithm>
#include <vector>

NAMESPACE_BEGIN(mitsuba)
/**
* Merge multiple lists of nodes into one. Remove duplicates.
*/
template <typename Float>
std::vector<Float> merge_nodes(const std::vector<std::vector<Float>> &nodes) {
    if (nodes.empty())
        return {};
    if (nodes.size() == 1)
        return nodes[0];

    std::vector<Float> merged_nodes;
    for (const auto &n : nodes)
        merged_nodes.insert(merged_nodes.end(), n.begin(), n.end());
    std::sort(merged_nodes.begin(), merged_nodes.end());
    // std::unique moves indeterminate duplicates at the end of the range
    // which are erased by std::vector::erase.
    merged_nodes.erase(std::unique(merged_nodes.begin(), merged_nodes.end()),
                       merged_nodes.end());

    return std::move(merged_nodes);
}

NAMESPACE_END(mitsuba)
