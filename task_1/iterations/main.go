package main

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

type Task struct {
	labio.Parameters
	ID         string        `json:"id"`
	Matrix     linalg.Matrix `json:"matrix"`
	RHS        []float64     `json:"rhs"`
	Comparison []float64     `json:"comparison_epsilons"`
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
		f, err := linalg.FactorLU(t.Matrix, false)
		if err != nil {
			return err
		}
		reference := linalg.SolveLU(f, t.RHS, false)
		fmt.Printf("task %s\nLU reference = %.12g\nmethod       epsilon      iterations   reference error   residual\n", t.ID, reference)
		counts := []int{}
		for _, epsilon := range append([]float64{t.Epsilon}, t.Comparison...) {
			for _, seidel := range []bool{false, true} {
				name := "jacobi"
				if seidel {
					name = "seidel"
				}
				r, err := linalg.LinearIterations(t.Matrix, t.RHS, epsilon, t.MaxIterations, seidel)
				if err != nil {
					return fmt.Errorf("%s/%s: %w", t.ID, name, err)
				}
				last := r.History[len(r.History)-1]
				error := linalg.Distance(r.X, reference)
				fmt.Printf("%-10s %10.3g %10d %17.5g %12.5g\n", name, epsilon, last.Iteration, error, last.Residual)
				if epsilon == t.Epsilon {
					fmt.Printf("  x = %.12g; certified bound %.3g, ||T^%d||inf = %.6g\n", r.X, last.Estimate, r.Block, r.Contraction)
					counts = append(counts, last.Iteration)
					if o.ExtraInfo {
						labio.PrintHistory(name, r.History, reference)
					}
				}
			}
		}
		fmt.Printf("at epsilon=%g: Jacobi %d iterations, Seidel %d iterations\n", t.Epsilon, counts[0], counts[1])

	}
	return nil
}
func main() { labio.Main(run) }
