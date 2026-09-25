# Re-resolve every Julia environment that path-devs the Epicycle subpackages.
# Run after adding/removing a dep in any subpackage `Project.toml`.
#
#   julia --startup-file=no ci/resolve_all.jl

using Pkg

const ENVS = (
    raw"c:\Users\steve\Dev\epicycle-dev",
    raw"c:\Users\steve\Dev\Epicycle",
)

for env in ENVS
    println("--- resolving $env ---")
    Pkg.activate(env)
    Pkg.resolve()
end

println("resolve_all: done")
