# Contributing

## Code Formatting

This project uses [Runic.jl](https://github.com/fredrikekre/Runic.jl) for consistent code formatting across all Julia files.

### Installing Runic.jl

Since Runic.jl is not available in the General registry, it needs to be installed in a separate environment:

```shell
julia --project=@runic -e 'using Pkg; Pkg.add(url = "https://github.com/fredrikekre/Runic.jl")'
```

### Formatting Code

To format all Julia files in the project:

```shell
julia --project=@runic -e 'using Runic; exit(Runic.main(["--inplace", "."]))'
```

This will recursively format all `.jl` files in the project directory using Runic's fixed formatting rules.

### Optional: Shell Alias

For convenience, you can add this alias to your shell configuration (`.bashrc`, `.zshrc`, etc.):

```shell
alias runic="julia --project=@runic -e 'using Runic; exit(Runic.main(ARGS))'"
```

Then you can simply run:

```shell
runic --inplace .
```
