package main

import (
	"fmt"
	"math"
	"nm/internal/linalg"
)

func iteration(t Task) ([]linalg.Step, error) {
	if t.Q < 0 || t.Q >= 1 {
		return nil, fmt.Errorf("simple iteration requires 0 <= q < 1")
	}
	x := t.Initial
	var history []linalg.Step
	for k := 1; k <= t.MaxIterations; k++ {
		next := phi(x)
		if math.IsNaN(next) || next < t.Lower || next > t.Upper {
			return nil, fmt.Errorf("iteration %d left the contraction interval", k)
		}
		bound := t.Q / (1 - t.Q) * math.Abs(next-x)
		history = append(history, linalg.Step{Iteration: k, X: []float64{next}, Estimate: bound, Residual: math.Abs(f(next))})
		if bound <= t.Epsilon {
			return history, nil
		}
		x = next
	}
	return nil, fmt.Errorf("simple iteration did not converge in %d iterations", t.MaxIterations)
}

func newton(t Task) ([]linalg.Step, error) {
	x := t.Initial
	var history []linalg.Step
	for k := 1; k <= t.MaxIterations; k++ {
		d := derivative(x)
		if d == 0 {
			return nil, fmt.Errorf("zero derivative at iteration %d", k)
		}
		next := x - f(x)/d
		if math.IsNaN(next) || math.IsInf(next, 0) {
			return nil, fmt.Errorf("non-finite Newton approximation at iteration %d", k)
		}
		step := math.Abs(next - x)
		residual := math.Abs(f(next))
		history = append(history, linalg.Step{Iteration: k, X: []float64{next}, Estimate: step, Residual: residual})
		if step <= t.Epsilon && residual <= t.Epsilon {
			return history, nil
		}
		x = next
	}
	return nil, fmt.Errorf("Newton did not converge in %d iterations", t.MaxIterations)
}
