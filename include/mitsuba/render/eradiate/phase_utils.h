#pragma once

#include <mitsuba/core/fwd.h>
#include <drjit/array.h>
#include <algorithm>
#include <vector>

NAMESPACE_BEGIN(mitsuba)

/**
* Merge multiple node buffers into one. Remove duplicates.
*
* The concatenation/sort/unique is performed on the host: JIT buffers are
* evaluated and copied to a temporary scalar array before being merged and
* loaded back into a single \c FloatStorage.
*/
template <typename FloatStorage>
FloatStorage merge_envelope_nodes(const std::vector<FloatStorage> &lists) {
    using ScalarFloat = dr::scalar_t<FloatStorage>;

    if (lists.empty())
        return FloatStorage();
    if (lists.size() == 1)
        return lists[0];

    std::vector<ScalarFloat> merged;
    for (const auto &l : lists) {
        size_t w   = dr::width(l);
        size_t off = merged.size();
        merged.resize(off + w);
        dr::store(merged.data() + off, l);
    }

    std::sort(merged.begin(), merged.end());
    // std::unique moves indeterminate duplicates at the end of the range
    // which are erased by std::vector::erase.
    merged.erase(std::unique(merged.begin(), merged.end()), merged.end());

    return dr::load<FloatStorage>(merged.data(), merged.size());
}

NAMESPACE_END(mitsuba)
