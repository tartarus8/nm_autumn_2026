package main

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

type Task struct {
	labio.Parameters
	ID     string        `json:"id"`
	Matrix linalg.Matrix `json:"matrix"`
}

func run(o labio.Options) error {
	tasks, err := labio.ReadTasks[Task](o.Input)
	if err != nil {
		return err
	}
	for _, t := range tasks {
		if err = t.Check(); err != nil {
			return err
		}
		fmt.Printf("task %s; epsilon=%g\n", t.ID, t.Epsilon)
		q, r := linalg.HouseholderQR(t.Matrix)
		linalg.PrintMatrix("A", t.Matrix)
		linalg.PrintMatrix("Q", q)
		linalg.PrintMatrix("R", r)
		product := linalg.Multiply(q, r)
		result, err := linalg.QREigenvalues(t.Matrix, t.Epsilon, t.MaxIterations)
		if err != nil {
			return err
		}
		fmt.Printf("eigenvalues = %.12g\nQR iterations = %d\n", result.Values, len(result.History))
		linalg.PrintMatrix("Q*R", product)
		fmt.Printf("max |QR-A| = %.3e; max |Q^T Q-I| = %.3e\n", linalg.MaxDifference(product, t.Matrix), linalg.MaxDifference(linalg.Multiply(linalg.Transpose(q), q), linalg.Identity(len(q))))
		if o.ExtraInfo {
			for _, step := range result.History {
				fmt.Printf("iteration %d: block size %d, %s shift, last subdiagonal %.5g\n", step.Iteration, step.Size, step.Shift, step.Subdiagonal)
			}
		}
		fmt.Println()
	}
	return nil
}
func main() { labio.Main(run) }
