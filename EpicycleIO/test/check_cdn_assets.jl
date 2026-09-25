# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Do the third-party browser libraries the dashboard depends on still resolve?
#
#   julia --project=<env> EpicycleIO/test/check_cdn_assets.jl
#
# EpicycleIO draws with Plotly and Cesium, and neither is vendored or installed. `dashboard.html`
# pulls both from jsdelivr at page-load time, pinned to exact versions. Nothing else in the test
# suite touches a browser: the 335 assertions exercise JSON and CZML generation, so if a pinned
# version were withdrawn from the CDN every user's plots and trajectories would go blank and the
# suite would stay green.
#
# This is deliberately NOT included from runtests.jl. It needs network, and a registered package
# whose tests require network fails for anyone testing offline. CI runs it as its own step, and
# a release should run it by hand.
#
# The URLs are read out of dashboard.html rather than repeated here, so the check cannot drift
# from the page it is checking.

using Downloads

const DASHBOARD = joinpath(@__DIR__, "..", "src", "assets", "dashboard.html")

"""Every absolute http(s) URL the dashboard loads, in the order it loads them."""
function dashboard_urls(path::AbstractString)
    html = read(path, String)
    urls = String[]
    # `src="…"`, `href="…"`, and the base-plus-suffix form the Cesium loader uses.
    for m in eachmatch(r"(?:src|href)\s*=\s*\"(https?://[^\"]+)\"", html)
        push!(urls, m[1])
    end
    base = match(r"CESIUM_BASE\s*=\s*'([^']+)'", html)
    if base !== nothing
        for m in eachmatch(r"CESIUM_BASE\s*\+\s*'([^']+)'", html)
            push!(urls, base[1] * m[1])
        end
    end
    return unique(urls)
end

function main()
    urls = dashboard_urls(DASHBOARD)
    isempty(urls) && (println("No URLs found in dashboard.html — has it been restructured?"); return 2)

    println("Third-party browser assets, from ", relpath(DASHBOARD), "\n")
    bad = String[]
    for url in urls
        ok, note = try
            # HEAD is enough and avoids pulling ~10 MB on every CI run.
            io = IOBuffer()
            r = Downloads.request(url; method = "HEAD", output = io, throw = false)
            (r.status == 200, string(r.status))
        catch e
            (false, sprint(showerror, e))
        end
        println("  ", ok ? "ok    " : "FAIL  ", rpad(note, 8), url)
        ok || push!(bad, url)
    end

    println()
    if isempty(bad)
        println(length(urls), " asset(s) reachable. A user with network gets a working dashboard.")
        return 0
    end
    println(length(bad), " asset(s) did NOT resolve. Every user's dashboard is affected, not just")
    println("this machine — the version is pinned in dashboard.html and has to be moved there.")
    return 1
end

exit(main())
