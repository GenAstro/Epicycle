# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Leap-second lists for tests, in the published format: every change through 2017, the update and
# expiry lines, and the #h hash computed as the published file computes it. The hash rule is held
# against a real list in test_correctness_leap_seconds.jl.

using SHA: sha1
using AstroEpochs: _BUILTIN_LEAP_SECONDS, cal2jd

const _LIST_NTP_EPOCH_JD = 2415020.5                 # 1900-01-01, the list's time origin

_list_ntp(jd) = round(Int, (jd - _LIST_NTP_EPOCH_JD) * 86400)

function _list_builtin_entries()
    return [(_list_ntp(sum(cal2jd(y, m, d))), Int(delta)) for (y, m, d, delta) in _BUILTIN_LEAP_SECONDS]
end

_list_hash(updated, expires, entries) =
    bytes2hex(sha1(string(updated) * string(expires) * join(string(a) * string(b) for (a, b) in entries)))

"""
    _valid_leap_list(expires_ntp; extra = [], entries = nothing, file_hash = nothing)

The text of a list. `extra` adds entries after 2017; `entries` replaces them all. The hash is computed
from what is written, so a test can make a list that hashes correctly and is still wrong.
"""
function _valid_leap_list(expires_ntp; extra = Tuple{Int, Any}[], entries = nothing, file_hash = nothing)
    updated = 3992312697
    rows = entries === nothing ? [_list_builtin_entries(); extra] : entries
    digest = file_hash === nothing ? _list_hash(updated, expires_ntp, rows) : file_hash
    grouped = join((digest[i:min(i + 7, end)] for i in 1:8:length(digest)), " ")
    return string("#\tA test list in the published format\n",
                  "#\$\t", updated, "\n",
                  "#@\t", expires_ntp, "\n",
                  join((string(a, "\t", b, "\t# entry\n") for (a, b) in rows)),
                  "#h\t", grouped, "\n")
end
