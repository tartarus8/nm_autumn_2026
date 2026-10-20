package nonlinear

import (
	"fmt"
	"math"
	"nm/internal/linalg"
)

type Problem struct {
	F                     Vector
	Jacobian              Jacobian
	Phi                   Vector
	PhiJacobian           Jacobian
	Initial, Lower, Upper []float64
	Q                     float64
}

func finite(x []float64) bool {
	for _, v := range x {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return false
		}
	}
	return true
}

func inside(x, lower, upper []float64) bool {
	for i, v := range x {
		if v < lower[i] || v > upper[i] {
			return false
		}
	}
	return true
}

func positive(x []float64) error {
	for _, v := range x {
		if v <= 0 {
			return fmt.Errorf("converged to a non-positive solution; choose another starting region")
		}
	}
	return nil
}

func Solve(p Problem, epsilon float64, limit int, newton bool) ([]linalg.Step, error) {
	if !newton && !(p.Q >= 0 && p.Q < 1) {
		return nil, fmt.Errorf("fixed-point iteration requires a supplied contraction bound 0 <= q < 1")
	}
	x := append([]float64(nil), p.Initial...)
	if !inside(x, p.Lower, p.Upper) {
		return nil, fmt.Errorf("initial approximation is outside the configured region")
	}
	fx := p.F(x)
	if !finite(fx) {
		return nil, fmt.Errorf("non-finite function at initial approximation")
	}
	history := []linalg.Step{{Iteration: 0, X: x, Estimate: math.Inf(1), Residual: linalg.Norm(fx)}}
	if linalg.Norm(fx) == 0 {
		history[0].Estimate = 0
		return history, positive(x)
	}
	for k := 1; k <= limit; k++ {
		var next []float64
		if newton {
			jac := p.Jacobian(x)
			for _, row := range jac {
				if !finite(row) {
					return history, fmt.Errorf("non-finite Jacobian at iteration %d", k)
				}
			}
			f, err := linalg.FactorLU(jac, false)
			if err != nil {
				return history, fmt.Errorf("Newton iteration %d: %w", k, err)
			}
			delta := linalg.SolveLU(f, fx, false)
			next = make([]float64, len(x))
			for i := range x {
				next[i] = x[i] - delta[i]
			}
		} else {
			jac := p.PhiJacobian(x)
			for _, row := range jac {
				if !finite(row) {
					return history, fmt.Errorf("non-finite iteration Jacobian at iteration %d", k)
				}
			}
			if norm := linalg.NormInf(jac); norm > p.Q+1e-12 {
				return history, fmt.Errorf("iteration %d: ||Phi'(x)||inf=%g exceeds supplied q=%g", k, norm, p.Q)
			}
			next = p.Phi(x)
		}
		if !finite(next) {
			return history, fmt.Errorf("non-finite approximation at iteration %d", k)
		}
		if !inside(next, p.Lower, p.Upper) {
			return history, fmt.Errorf("iteration %d left the configured region; review the initial guess, region, or iteration mapping", k)
		}
		values := p.F(next)
		if !finite(values) {
			return history, fmt.Errorf("non-finite function at iteration %d", k)
		}
		step := linalg.Distance(next, x)
		estimate := step
		if !newton {
			estimate = p.Q / (1 - p.Q) * step
		}
		residual := linalg.Norm(values)
		history = append(history, linalg.Step{Iteration: k, X: next, Estimate: estimate, Residual: residual})
		if estimate <= epsilon && (!newton || residual <= epsilon) {
			return history, positive(next)
		}
		x, fx = next, values
	}
	return history, fmt.Errorf("nonlinear solver did not converge within %d iterations", limit)
}
