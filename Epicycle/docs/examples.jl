# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Example pages, generated from the example scripts.
#
# An example has one copy: the script in `examples/`, which ships with the package and runs on its
# own. Literate turns that script into a documentation page, so the page cannot drift from the code
# and the nav shows what the package can do.
#
# A script marks the prose that belongs on the page with `#'`, the way a docstring marks itself.
# Everything else stays what it is: an ordinary `#` comment labelling a step, which the page shows
# inside the code block. The licence header is dropped from the page.
#
# Running the examples is a separate decision from rendering them, because they take minutes and a
# documentation change does not touch them:
#
#     julia --project=<env> docs/make.jl                          render only, no output shown
#     EPICYCLE_DOC_EXAMPLES=all julia --project=<env> docs/make.jl     run every example
#     EPICYCLE_DOC_EXAMPLES=Ex_OrbitDetermination,Ex_SimpleTarget ...  run those two
#
# What an example prints is verified by the regression harness, which checks the numbers rather
# than that the script completed. The switch here is for a page that shows its output.

using Literate

"Which examples this build runs: `all`, `none`, or the script names to run."
const RUN_EXAMPLES = get(ENV, "EPICYCLE_DOC_EXAMPLES", "none")

runs_example(name) = RUN_EXAMPLES == "all" ||
                     name in strip.(split(RUN_EXAMPLES, ","; keepempty = false))

# `#'` opens prose for the page. An ordinary comment becomes `##`, which Literate keeps inside the
# code block rather than promoting to prose, and the licence header goes entirely.
function literate_source(str)
    out = IOBuffer()
    for line in split(str, '\n')
        stripped = lstrip(line)
        if occursin(r"^#\s*(Copyright \(C\)|SPDX-License-Identifier)", stripped)
            continue
        elseif startswith(stripped, "#' ")
            println(out, replace(line, "#' " => "# "; count = 1))
        elseif startswith(stripped, "#'")
            println(out, replace(line, "#'" => "#"; count = 1))
        elseif startswith(stripped, "#")
            println(out, replace(line, "#" => "##"; count = 1))
        else
            println(out, line)
        end
    end
    return String(take!(out))
end

"""
    example_pages(groups, examples_dir, out_dir) -> Vector

Generate a page for each script named in `groups` and return the `pages` entry for `makedocs`.

`groups` is a vector of `"Group title" => ["Page title" => "script_name", …]`, where the script is
`examples_dir/<script_name>.jl`. Pages are written to `out_dir`, which is emptied first so a
renamed script leaves nothing behind.
"""
function example_pages(groups, examples_dir, out_dir)
    isdir(out_dir) && rm(out_dir; recursive = true)
    mkpath(out_dir)

    nav = Any[]
    for (group, entries) in groups
        pages = Any[]
        for (title, script) in entries
            path = joinpath(examples_dir, script * ".jl")
            isfile(path) || error("example_pages: no script at $path")
            execute = runs_example(script)
            @info "example page: $script" execute
            Literate.markdown(path, out_dir;
                              name      = script,
                              execute   = execute,
                              flavor    = Literate.DocumenterFlavor(),
                              codefence = "```julia" => "```",   # rendered, never re-run
                              preprocess = literate_source,
                              credit    = false)
            push!(pages, title => joinpath(basename(out_dir), script * ".md"))
        end
        push!(nav, group => pages)
    end
    return nav
end
