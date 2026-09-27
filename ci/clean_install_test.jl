# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Install Epicycle the way a user does, and run an example.
#
# This tests the published registries, not the working tree, so its result is independent of the
# commit that triggered it. Nothing is cached: the point is a cold install, including the SPICE
# kernel download, on a machine that has never seen Epicycle.
#
#     julia --project=@epicycle-install ci/clean_install_test.jl
#
# Throwing is how it fails, so the exit code carries the result.

using Pkg

println("\n=== julia ===")
println(VERSION, "  ", Sys.MACHINE)
println(DEPOT_PATH[1])

println("\n=== registries ===")
# Both registries in one call. Julia fetches General on its own only while no registry is
# present at all, so adding ours first would leave General out and every dependency
# unresolvable. This is the functional form of `registry add General <url>`; Pkg warns
# against the REPL mode from a script.
Pkg.Registry.add([
    RegistrySpec(name = "General"),
    RegistrySpec(url  = "https://github.com/GenAstro/GenAstroRegistry.git"),
])
println("registries: ", [r.name for r in Pkg.Registry.reachable_registries()])

println("\n=== install ===")
Pkg.add("Epicycle")

println("\n=== versions resolved ===")
Pkg.status()

# A path instead of a version means something resolved to a checkout rather than the registry,
# which would make the whole run meaningless.
deps = Pkg.dependencies()
developed = [d.name for (_, d) in deps if d.is_tracking_path]
isempty(developed) || error("resolved to a local path, not the registry: ", join(developed, ", "))

println("\n=== load, and download the kernels ===")
using Epicycle
println("Epicycle ", pkgversion(Epicycle))

println("\n=== run an example ===")
# The viewer warns rather than throws on a machine with no browser, so this is headless-safe.
Epicycle.run_example("Ex_GettingStarted"; echo = false)


println("\n=== every documented install path ===")
# Each package README tells a user to add that package on its own. Install each one that way, in
# its own temporary environment, and check it resolved to the newest version the registry offers
# rather than to an older copy General still carries. The expected version is read from the
# registry, so this needs no editing when something is released.
#
# Resolution is what can fail here; compilation already succeeded for all of these as part of the
# Epicycle install above. Turning auto-precompilation off keeps each check to a resolve and a
# manifest write rather than a rebuild in every temporary environment.
withenv("JULIA_PKG_PRECOMPILE_AUTO" => "0") do

gen = only(filter(r -> r.name == "GenAstroRegistry", Pkg.Registry.reachable_registries()))
for (_, entry) in sort(collect(gen.pkgs); by = x -> x[2].name)
    want = maximum(keys(Pkg.Registry.registry_info(entry).version_info))
    Pkg.activate(temp = true)
    Pkg.add(entry.name)
    got = first(d.version for (_, d) in Pkg.dependencies() if d.name == entry.name)
    got == want || error("$(entry.name) resolved to $got; the registry offers $want")
    println("  ", rpad(entry.name, 16), got, "  OK")
end

end

println("\n>>> CLEAN INSTALL OK")
