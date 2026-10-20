package main

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

func run(extraInfo bool) error {
	for _, task := range tasks {
		fmt.Printf("task %d\n", task.ID)
		if extraInfo {
			linalg.PrintMatrix("A", task.Matrix)
		}
		f, err := linalg.FactorLU(task.Matrix, extraInfo)
		if err != nil {
			return fmt.Errorf("task %d: %w", task.ID, err)
		}
		x := linalg.SolveLU(f, task.RHS, extraInfo)
		det := linalg.DeterminantLU(f)
		inverse := linalg.InverseLU(f, extraInfo)

		if extraInfo {
			linalg.PrintMatrix("L", f.L)
			linalg.PrintMatrix("U", f.U)
			linalg.PrintMatrix("P", f.P)
		}
		fmt.Printf("\nx = %.10g\ndet(A) = %.10g\n\n", x, det)
		linalg.PrintMatrix("A^-1", inverse)

		product := linalg.Multiply(task.Matrix, inverse)
		linalg.PrintMatrix("A * A^-1", product)
	}
	return nil
}

func main() {
	labio.Main(run)
}
