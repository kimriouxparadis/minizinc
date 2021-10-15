using PackageCompiler
using Pkg
Pkg.resolve()
Pkg.instantiate()
create_app("Exec", "Compiled", force = true)