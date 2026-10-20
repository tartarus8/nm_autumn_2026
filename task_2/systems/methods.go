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
		for i, value := range next {
			if math.IsNaN(value) || value < t.Lower[i] || value > t.Upper[i] {
				return nil, fmt.Errorf("iteration %d left the contraction region", k)
			}
		}
		bound := t.Q / (1 - t.Q) * linalg.Distance(next, x)
		history = append(history, linalg.Step{Iteration: k, X: next, Estimate: bound, Residual: linalg.Norm(f(next))})
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
		factors, err := linalg.FactorLU(jacobian(x), false)
		if err != nil {
			return nil, err
		}
		delta := linalg.SolveLU(factors, f(x), false)
		next := make([]float64, len(x))
		for i := range x {
			next[i] = x[i] - delta[i]
			if math.IsNaN(next[i]) || math.IsInf(next[i], 0) {
				return nil, fmt.Errorf("non-finite Newton approximation at iteration %d", k)
			}
		}
		step := linalg.Distance(next, x)
		residual := linalg.Norm(f(next))
		history = append(history, linalg.Step{Iteration: k, X: next, Estimate: step, Residual: residual})
		if step <= t.Epsilon && residual <= t.Epsilon {
			return history, nil
		}
		x = next
	}
	return nil, fmt.Errorf("Newton did not converge in %d iterations", t.MaxIterations)
}
