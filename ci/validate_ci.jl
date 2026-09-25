# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

# Pre-flight checks for the CI configuration. Seconds, not a full cycle.
#
#     julia --startup-file=no --project=. ci/validate_ci.jl
#
# `ci/test_locally.jl` builds an isolated depot and precompiles the whole stack from nothing, so it
# takes about forty-five minutes to tell you that a package name is missing from a list. Every
# check here is one that cost a full cycle on 2026-09-16, and every one of them is answerable by
# reading files.
#
# This loads no packages beyond stdlib and resolves nothing, so it is safe to run at any time and
# against any branch. It reports rather than fixes.

using TOML
using Pkg

const REPO_ROOT = dirname(@__DIR__)

const RED    = "❌"
const GREEN  = "✅"
const YELLOW = "⚠️ "

const problems = String[]
const warnings = String[]

fail(msg)  = (push!(problems, msg); println("  $RED $msg"))
warn_(msg) = (push!(warnings, msg); println("  $YELLOW $msg"))
pass(msg)  = println("  $GREEN $msg")

# ── Reading the hand-maintained lists ────────────────────────────────────────
#
# The CI scripts hold package names in plain array literals. Parsing the Julia would be exact and
# brittle; a regex over the literal is neither, but it fails loudly when a list is renamed, which
# is the behaviour that matters — a silently empty list would make every check below vacuous.

"""The workspace members, which is the one list the others are derived from."""
workspace_projects() =
    TOML.parsefile(joinpath(REPO_ROOT, "Project.toml"))["workspace"]["projects"]

"""
Extract the names assigned to `var` in `path`.

Two forms are in use. A literal array is read directly. A list derived from the workspace — either
`WORKSPACE_PACKAGES` or a direct read of `["workspace"]["projects"]` — resolves to the workspace
members, with `Epicycle` dropped where the expression filters it out. Reading only the literal form
reports every derived list as missing, and a list that cannot be read looks the same as an empty
one, so every downstream check then reports that nothing is built or tested.
"""
function list_names(path::AbstractString, var::AbstractString)
    isfile(path) || return nothing
    src = read(path, String)
    m = match(Regex(var * raw"\s*=\s*\[(.*?)\]", "s"), src)
    if m !== nothing
        return [String(x.captures[1]) for x in eachmatch(r"[:\"]([A-Za-z][A-Za-z0-9_]*)\"?", m.captures[1])]
    end
    d = match(Regex(var * raw"\s*=(.{0,400})", "s"), src)
    d === nothing && return nothing
    rhs = d.captures[1]
    occursin("WORKSPACE_PACKAGES", rhs) || occursin("\"workspace\"", rhs) || return nothing
    names = workspace_projects()
    occursin("!=(\"Epicycle\")", rhs) && (names = filter(!=("Epicycle"), names))
    return String.(names)
end

# (label, file, variable, what the list means)
const LISTS = [
    ("setup:develop",  "ci/setup_environment.jl",            "packages",             "developed into the project"),
    ("build:umbrella", "ci/build_epicycle.jl",               "packages_to_check",    "must be re-exported by Epicycle"),
    ("build:docs",     "ci/build_epicycle.jl",               "packages_to_document", "documentation is built"),
    ("test:packages",  "ci/test_epicycle.jl",                "packages",             "tests are run"),
    ("test:covscan",   "ci/test_epicycle.jl",                "packages_to_check",    "scanned for .cov files"),
    ("coverage",       "ci/generate_coverage.jl",            "packages",             "coverage is collected"),
    ("suite",          "Epicycle/util/test_all_packages.jl", "packages",             "suite runner"),
]

println("\n", "="^72)
println("CI configuration pre-flight")
println("="^72)

const lists = Dict{String,Vector{String}}()
println("\n▸ Reading the package lists")
for (label, file, var, _) in LISTS
    names = list_names(joinpath(REPO_ROOT, file), var)
    if names === nothing
        fail("could not read `$var` in $file — renamed or restructured?")
        continue
    end
    lists[label] = names
    pass("$label: $(length(names)) packages  ($file → $var)")
end

is_package(p) = isfile(joinpath(REPO_ROOT, p, "Project.toml")) &&
                isfile(joinpath(REPO_ROOT, p, "src", "$p.jl"))

# ── 1. Every listed package exists on disk ───────────────────────────────────

println("\n▸ Every listed package exists")
for (label, names) in sort(collect(lists); by = first)
    absent = filter(!is_package, names)
    isempty(absent) ? pass("$label — all present") :
                      fail("$label names packages not on disk: $(join(absent, ", "))")
end

# ── 2. The umbrella load check matches what the umbrella re-exports ──────────
#
# `packages_to_check` asks `isdefined(Main, :X)` after `using Epicycle`, so it can only hold names
# the umbrella actually re-exports. Adding a package the umbrella does not depend on makes CI fail
# with "X failed to load" however healthy X is. That cost a cycle on 2026-09-16.

println("\n▸ Umbrella load check matches the umbrella's re-exports")
let umbrella = joinpath(REPO_ROOT, "Epicycle", "src", "Epicycle.jl")
    if isfile(umbrella)
        reexported = Set(String(m.captures[1]) for m in
                         eachmatch(r"@reexport\s+using\s+([A-Za-z][A-Za-z0-9_]*)", read(umbrella, String)))
        listed = get(lists, "build:umbrella", String[])
        bad = filter(n -> !(n in reexported), listed)
        isempty(bad) ? pass("all $(length(listed)) names are re-exported by Epicycle") :
            fail("`packages_to_check` names $(join(bad, ", ")), which Epicycle does not @reexport — " *
                 "isdefined(Main, ...) will be false and CI will report a false failure")
    else
        warn_("no Epicycle/src/Epicycle.jl — skipped")
    end
end

# ── 3. Root [deps] entries can survive a fresh resolve ───────────────────────
#
# A fresh CI checkout has no manifest, so `Pkg.instantiate()` resolves every [deps] entry from the
# registry. A package that lives in this repo and is not yet registered cannot be a plain [deps]
# entry — it has to arrive through `Pkg.develop`, which writes a path into the manifest. Committing
# such an entry breaks CI at the first step. That cost a cycle on 2026-09-16.
#
# Ask Pkg rather than reading the registry directory: a registry is normally stored compressed, so
# looking for unpacked Package.toml files reports "cannot tell" on a perfectly healthy machine.

println("\n▸ Root [deps] entries survive a fresh resolve")
const REGISTERED = try
    Set(pkg.name for r in Pkg.Registry.reachable_registries() for (_, pkg) in r.pkgs)
catch
    nothing
end

let root_project = TOML.parsefile(joinpath(REPO_ROOT, "Project.toml")),
    root_deps    = sort(collect(keys(get(root_project, "deps", Dict{String,Any}()))))

    global ROOT_DEPS = root_deps
    local local_pkgs = filter(is_package, root_deps)
    if REGISTERED === nothing
        warn_("no reachable registry — cannot check whether repo-local deps are registered")
    elseif isempty(local_pkgs)
        pass("no repo-local packages are plain [deps] entries")
    else
        # An unregistered repo-local dep is fine when setup develops the workspace before it
        # instantiates, because develop writes a path into the manifest and the resolve then has
        # somewhere to look. The invariant is the ordering, so check that rather than the symptom.
        local setup_src   = read(joinpath(REPO_ROOT, "ci", "setup_environment.jl"), String)
        local dev_at      = findfirst("Pkg.develop", setup_src)
        local inst_at     = findfirst("Pkg.instantiate", setup_src)
        local develops_first = dev_at !== nothing &&
                               (inst_at === nothing || first(dev_at) < first(inst_at))
        local workspace = workspace_projects()
        for dep in local_pkgs
            if dep in REGISTERED
                pass("$dep is registered — a plain [deps] entry is fine")
            elseif dep in workspace && develops_first
                pass("$dep is unregistered but is a workspace member, and setup develops before " *
                     "it instantiates — the resolve finds it by path")
            elseif dep in workspace
                fail("$dep is in root [deps], is NOT registered, and ci/setup_environment.jl " *
                     "instantiates before it develops — the resolve looks $dep up in the " *
                     "registries and fails with \"expected package $dep to be registered\". " *
                     "Move the develop loop above Pkg.instantiate()")
            else
                fail("$dep is in root [deps] but is NOT registered and is not a workspace " *
                     "member — a fresh Pkg.instantiate() fails with \"expected package $dep to " *
                     "be registered\"; it must arrive via Pkg.develop until it is registered")
            end
        end
    end
end

# ── 4. The manifest, if present, holds every direct dependency ───────────────
#
# `setup_environment.jl` calls instantiate before resolve, so a manifest missing a direct dep kills
# the run before the resolve that would have fixed it. `Revise` had been missing since 2026-09-13.

println("\n▸ Manifest is consistent with Project.toml")
let manifest = joinpath(REPO_ROOT, "Manifest.toml")
    if !isfile(manifest)
        pass("no manifest — instantiate resolves from scratch, as CI does")
    else
        entries = Set(keys(get(TOML.parsefile(manifest), "deps", Dict{String,Any}())))
        absent  = filter(d -> !(d in entries), ROOT_DEPS)
        isempty(absent) ? pass("every direct dependency is in the manifest") :
            fail("manifest is missing direct deps: $(join(sort(absent), ", ")) — " *
                 "Pkg.instantiate() refuses before Pkg.resolve() runs; delete Manifest.toml or resolve first")
    end
end

# ── 5. Executed documentation blocks can find what they import ───────────────
#
# Documenter runs `@example`, `@repl` and `jldoctest` blocks in the *active* project, which in CI is
# the repository root. A block importing something the root project does not declare fails the docs
# build — and builds fine against a warm development environment that happens to have it, which is
# how it reaches CI unnoticed. ForwardDiff did exactly that on 2026-09-16.

println("\n▸ Executed doc blocks import only what the root project provides")
let available = Set(ROOT_DEPS) ∪ Set(readdir(Sys.STDLIB)),
    blocks = 0, files = Set{String}()

    for pkg in sort(get(lists, "build:docs", String[]))
        srcdir = joinpath(REPO_ROOT, pkg, "docs", "src")
        isdir(srcdir) || continue
        push!(available, pkg)                       # a package's own docs may import itself
        for (root, _, fs) in walkdir(srcdir), f in fs
            endswith(f, ".md") || continue
            text = read(joinpath(root, f), String)
            # Only blocks Documenter executes. A plain ```julia fence is never run.
            for blk in eachmatch(r"```(@example|@repl|jldoctest)[^\n]*\n(.*?)```"s, text)
                blocks += 1
                push!(files, relpath(joinpath(root, f), REPO_ROOT))
                for line in eachmatch(r"^[ \t]*using[ \t]+([A-Za-z][\w, ]*)"m, blk.captures[2])
                    for raw in split(line.captures[1], ",")
                        name = strip(raw)
                        (isempty(name) || name in available) && continue
                        fail("$pkg docs import `$name`, which the root project does not provide — " *
                             relpath(joinpath(root, f), REPO_ROOT))
                    end
                end
            end
        end
    end
    blocks == 0 ? warn_("no executed doc blocks found — are all fences plain ```julia?") :
                  pass("scanned $blocks executed block(s) across $(length(files)) file(s)")
end

# ── 5b. The scripts CI runs import only what the root project provides ───────
#
# CI activates the repository root and then includes these scripts, so a package they import has
# to be in the root [deps]. The warm `epicycle-dev` environment has more in it than the root does,
# which is why one of these builds locally and fails in CI. `Literate` was used by
# `Epicycle/docs/examples.jl` and declared nowhere, and the check above did not see it because it
# reads executed doc blocks and this is an ordinary `using` in a script.

println("
▸ Scripts CI runs import only what the root project provides")
let available = Set(ROOT_DEPS) ∪ Set(readdir(Sys.STDLIB)) ∪ Set(workspace_projects()),
    scanned = 0, bad = 0

    scripts = String[]
    append!(scripts, [joinpath("ci", f) for f in readdir(joinpath(REPO_ROOT, "ci")) if endswith(f, ".jl")])
    for pkg in sort(workspace_projects())
        d = joinpath(REPO_ROOT, pkg, "docs")
        isdir(d) || continue
        append!(scripts, [joinpath(pkg, "docs", f) for f in readdir(d) if endswith(f, ".jl")])
    end

    for rel in scripts
        scanned += 1
        for line in eachline(joinpath(REPO_ROOT, rel))
            startswith(lstrip(line), "#") && continue
            m = match(r"^[ 	]*(?:using|import)[ 	]+([A-Za-z][\w.:, ]*)", line)
            m === nothing && continue
            for raw in split(m.captures[1], ",")
                name = String(first(split(strip(raw), ('.', ':', ' '); keepempty = false), 1)[1])
                (isempty(name) || name in available) && continue
                fail("$rel imports `$name`, which the root project does not provide — " *
                     "CI activates the repository root, so add it to the root Project.toml")
                bad += 1
            end
        end
    end
    bad == 0 && pass("scanned $scanned script(s), every import is available at the root")
end

# ── 6. Nothing is silently skipped ───────────────────────────────────────────

println("\n▸ Nothing is silently skipped")
let skipped = 0
    for pkg in sort(filter(is_package, readdir(REPO_ROOT)))
        if isfile(joinpath(REPO_ROOT, pkg, "docs", "make.jl")) && !(pkg in get(lists, "build:docs", String[]))
            warn_("$pkg has docs/make.jl but is not in packages_to_document — its docs never build")
            skipped += 1
        end
        if isdir(joinpath(REPO_ROOT, pkg, "test")) && !(pkg in get(lists, "suite", String[]))
            warn_("$pkg has a test/ directory but is not in the suite runner — its tests never run")
            skipped += 1
        end
    end
    skipped == 0 && pass("every package on disk is documented and tested by CI")
end

# ── Verdict ──────────────────────────────────────────────────────────────────

println("\n", "="^72)
if isempty(problems) && isempty(warnings)
    println("$GREEN CI configuration looks sound. Safe to run ci/test_locally.jl.")
elseif isempty(problems)
    println("$GREEN No blocking problems. $(length(warnings)) warning(s) above.")
else
    println("$RED $(length(problems)) problem(s) would fail CI:")
    for p in problems
        println("   • ", first(split(p, " — ")))
    end
end
println("="^72)

exit(isempty(problems) ? 0 : 1)
