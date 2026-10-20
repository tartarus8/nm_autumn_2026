package main

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

type Task struct {
	ID       string    `json:"id"`
	Lower    []float64 `json:"lower"`
	Diagonal []float64 `json:"diagonal"`
	Upper    []float64 `json:"upper"`
	RHS      []float64 `json:"rhs"`
}

func run(o labio.Options) error {
	tasks, err := labio.ReadTasks[Task](o.Input)
	if err != nil {
		return err
	}
	for _, t := range tasks {
		x, p, q, err := linalg.ThreeDiag(t.Lower, t.Diagonal, t.Upper, t.RHS)
		if err != nil {
			return fmt.Errorf("%s: %w", t.ID, err)
		}
		fmt.Printf("task %s\np = %.10g\nq = %.10g\nx = %.10g\n", t.ID, p, q, x)
		if o.ExtraInfo {
			fmt.Printf("backward sweep: x[%d] = q[%d] = %.10g\n", len(x), len(x), x[len(x)-1])
			for i := len(x) - 2; i >= 0; i-- {
				fmt.Printf("x[%d] = (%.10g)*x[%d] + (%.10g) = %.10g\n", i+1, p[i], i+2, q[i], x[i])
			}
			fmt.Printf("max |Ax-b| = %.3e\n", linalg.Distance(linalg.ThreeDiagProduct(t.Lower, t.Diagonal, t.Upper, x), t.RHS))
		}
		fmt.Println()
	}
	return nil
}
func main() { labio.Main(run) }
