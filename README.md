# Numerical methods labs (Go)

Run from the repository root:

```sh
go run ./task_1/lu -extraInfo
go run ./task_1/ThreeDiag -extraInfo
go run ./task_1/iterations -extraInfo
go run ./task_1/rotations -extraInfo
go run ./task_1/qr -extraInfo

go run ./task_2/scalar -extraInfo
go run ./task_2/systems -extraInfo
```

- All required answers and verification products print by default. Add `-extraInfo` for intermediate steps.
- Each command reads `tests.json` beside its source. Keep the source directory available when running a compiled binary; do not use `-trimpath` with this default.
- All commands except `lu` also accept `-input FILE`.
- Results print to the terminal; `-extraInfo` adds steps and convergence tables. No result files are generated.
- Static Lab 2 starting-point pictures: [scalar](task_2/scalar/initial_guess.png) and [system](task_2/systems/initial_guess.png). These illustrate the bundled variant 1 inputs; replace them if you change the problems.
- Lab 2 uses provisional **variant 1** from the [course sheets](https://mainfo.ru/files/numeric-methods-lab-2.zip), not a confirmed personal assignment. The scalar task finds the smaller positive root; the system uses `a=2`. Formulas, initial values, reference solutions, and explanations of contraction bounds are in JSON.
- Lab 1 input files contain photo variants 5 and 6. QR uses a stricter tolerance of `1e-6`; the other iterative inputs use `0.01`. The additional required complex-pair example runs with `go run ./task_1/qr -input task_1/qr/examples/complex.json`.
- Shared algorithms are in `internal`; no external numerical libraries are required.

Run the tests:

```sh
go test ./...
```
