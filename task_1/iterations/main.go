package main

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

func run(extraInfo bool) error {
	for _, t := range tasks {
		if err := t.Check(); err != nil {
			return err
		}
		f, err := linalg.FactorLU(t.Matrix, false)
		if err != nil {
			return err
		}
		reference := linalg.SolveLU(f, t.RHS, false)
		fmt.Printf("task %s; epsilon=%g\nLU reference = %.12g\n", t.ID, t.Epsilon, reference)
		counts := []int{}
		for _, seidel := range []bool{false, true} {
			name := "jacobi"
			if seidel {
				name = "seidel"
			}
			r, err := linalg.LinearIterations(t.Matrix, t.RHS, t.Epsilon, t.MaxIterations, seidel)
			if err != nil {
				return fmt.Errorf("%s/%s: %w", t.ID, name, err)
			}
			last := r.History[len(r.History)-1]
			counts = append(counts, last.Iteration)
			fmt.Printf("\n%s: x = %.12g; iterations = %d\n", name, r.X, last.Iteration)
			fmt.Printf("reference error = %.3e; error bound = %.3e; residual = %.3e\n",
				linalg.Distance(r.X, reference), last.Estimate, last.Residual)
			if extraInfo {
				fmt.Printf("||T^%d||inf = %.6g\n\n", r.Block, r.Contraction)
				labio.PrintHistory(name, r.History, reference)
			}
		}
		fmt.Printf("\nat epsilon=%g: jacobi %d iterations, seidel %d iterations\n\n", t.Epsilon, counts[0], counts[1])
	}
	return nil
}

func main() {
	labio.Main(run)
}
