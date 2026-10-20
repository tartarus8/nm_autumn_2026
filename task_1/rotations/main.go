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
		r, err := linalg.Rotations(t.Matrix, t.Epsilon, t.MaxIterations)
		if err != nil {
			return err
		}
		fmt.Printf("task %s; epsilon=%g; rotations=%d\neigenvalues = %.12g\n", t.ID, t.Epsilon, len(r.History)-1, r.Values)
		lambda := linalg.Diagonal(r.Values)
		vt := linalg.Transpose(r.Vectors)
		for j, vector := range vt {
			fmt.Printf("v%d = %.12g\n", j+1, vector)
		}
		linalg.PrintMatrix("lambda", lambda)
		linalg.PrintMatrix("V (column i is the eigenvector for eigenvalue i)", r.Vectors)
		av, vl := linalg.Multiply(t.Matrix, r.Vectors), linalg.Multiply(r.Vectors, lambda)
		linalg.PrintMatrix("A*V", av)
		linalg.PrintMatrix("V*lambda", vl)
		fmt.Printf("max |AV-Vlambda| = %.3e; max |V^T V-I| = %.3e\n\n", linalg.MaxDifference(av, vl), linalg.MaxDifference(linalg.Multiply(vt, r.Vectors), linalg.Identity(len(t.Matrix))))
		if extraInfo {
			labio.PrintHistory("rotation (estimate: off-diagonal norm)", r.History, nil)
		}

	}
	return nil
}
func main() {
	labio.Main(run)
}
