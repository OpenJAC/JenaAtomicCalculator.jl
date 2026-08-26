# API: infrastructure

The modules JAC is built on: grids and orbitals, the nuclear model, the self-consistent field, the many-electron
basis, and the tabulated data behind them. They are published in full because a user writing a script reaches for
them directly -- `Nuclear.Model`, `Radial.Grid`, `AsfSettings` -- even though they carry no example files of their
own (Rule 16 of CLAUDE.md).

##  Atomic states and bases
```@autodocs ; canonical=false
Modules = [AtomicState]
Order   = [:type, :function]
```

##  Default settings and units
```@autodocs ; canonical=false
Modules = [Defaults]
Order   = [:type, :function]
```

##  Distributions
```@autodocs ; canonical=false
Modules = [Distribution]
Order   = [:type, :function]
```

##  Hamiltonian matrices
```@autodocs ; canonical=false
Modules = [Hamiltonian]
Order   = [:type, :function]
```

##  Many-electron bases and levels
```@autodocs ; canonical=false
Modules = [ManyElectron]
Order   = [:type, :function]
```

##  Nuclear models
```@autodocs ; canonical=false
Modules = [Nuclear]
Order   = [:type, :function]
```

##  Periodic table data
```@autodocs ; canonical=false
Modules = [PeriodicTable]
Order   = [:type, :function]
```

##  Radial grids and orbitals
```@autodocs ; canonical=false
Modules = [Radial]
Order   = [:type, :function]
```

##  Radial integrals
```@autodocs ; canonical=false
Modules = [RadialIntegrals]
Order   = [:type, :function]
```

##  Self-consistent field
```@autodocs ; canonical=false
Modules = [SelfConsistent]
Order   = [:type, :function]
```

##  Statistical tensors of an ensemble
```@autodocs ; canonical=false
Modules = [Statistical]
Order   = [:type, :function]
```
