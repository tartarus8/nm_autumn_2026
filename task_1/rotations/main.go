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
		r, err := linalg.Rotations(t.Matrix, t.Epsilon, t.MaxIterations)
		if err != nil {
			return err
		}
		fmt.Printf("task %s; epsilon=%g; rotations=%d\neigenvalues = %.12g\n", t.ID, t.Epsilon, len(r.History)-1, r.Values)
		lambda := linalg.Diagonal(r.Values)
		vt := linalg.Transpose(r.Vectors)
		for j, vector := range vt {
			fmt.Printf("v%d (lambda=%.12g) = %.12g\n", j+1, r.Values[j], vector)
		}
		linalg.PrintMatrix("Lambda", lambda)
		linalg.PrintMatrix("V (column j is the eigenvector for eigenvalue j)", r.Vectors)
		av, vl := linalg.Multiply(t.Matrix, r.Vectors), linalg.Multiply(r.Vectors, lambda)
		linalg.PrintMatrix("A*V", av)
		linalg.PrintMatrix("V*Lambda", vl)
		fmt.Printf("max |AV-VLambda| = %.3e; max |V^T V-I| = %.3e\n\n", linalg.MaxDifference(av, vl), linalg.MaxDifference(linalg.Multiply(vt, r.Vectors), linalg.Identity(len(t.Matrix))))
		if o.ExtraInfo {
			labio.PrintHistory("rotation (estimate: off-diagonal norm)", r.History, nil)
		}

	}
	return nil
}
func main() { labio.Main(run) }
