package linalg

import (
	"fmt"
	"math"
)

type IterativeResult struct {
	X           []float64
	History     []Step
	Contraction float64
	Block       int
}

// LinearIterations uses Jacobi, or Gauss-Seidel when seidel is true
func LinearIterations(a Matrix, b []float64, epsilon float64, limit int, seidel bool) (IterativeResult, error) {
	n := len(a)
	alpha, beta := Zeros(n, n), make([]float64, n)
	for i := range a {
		if a[i][i] == 0 {
			return IterativeResult{}, fmt.Errorf("zero diagonal at row %d: reorder equations first", i+1)
		}
		beta[i] = b[i] / a[i][i]
		for j := range a {
			if i != j {
				alpha[i][j] = -a[i][j] / a[i][i]
			}
		}
	}
	t := Clone(alpha)
	if seidel {
		t = Zeros(n, n)
		for i := range a {
			for j := range a {
				if j > i {
					t[i][j] = alpha[i][j]
				}
				for k := 0; k < i; k++ {
					t[i][j] += alpha[i][k] * t[k][j]
				}
			}
		}
	}
	power, block, q := Clone(t), 1, NormInf(t)
	for q >= 1 && block < limit && !math.IsInf(q, 0) {
		power = Multiply(power, t)
		block++
		q = NormInf(power)
	}
	if !(q < 1) {
		return IterativeResult{}, fmt.Errorf("could not find an iteration matrix within %d powers", limit)
	}
	x := append([]float64(nil), beta...)
	result := IterativeResult{Contraction: q, Block: block}
	result.History = append(result.History, Step{0, append([]float64(nil), x...), math.Inf(1), Distance(MatVec(a, x), b)})
	for iteration := 1; iteration <= limit; iteration++ {
		next := append([]float64(nil), beta...)
		for i := range a {
			for j := range a {
				value := x[j]
				if seidel && j < i {
					value = next[j]
				}
				next[i] += alpha[i][j] * value
			}
		}
		estimate := math.Inf(1)
		if iteration >= block {
			estimate = q / (1 - q) * Distance(next, result.History[iteration-block].X)
		}
		residual := Distance(MatVec(a, next), b)
		if math.IsNaN(residual) || math.IsInf(residual, 0) {
			return result, fmt.Errorf("non-finite iteration at step %d", iteration)
		}
		result.History = append(result.History, Step{iteration, next, estimate, residual})
		x = next
		if estimate <= epsilon {
			result.X = x
			return result, nil
		}
	}
	return result, fmt.Errorf("iteration did not converge within %d steps", limit)
}
