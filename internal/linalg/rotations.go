package linalg

import (
	"fmt"
	"math"
)

type EigenResult struct {
	Values  []float64
	Vectors Matrix
	History []Step
}

func Rotations(original Matrix, epsilon float64, limit int) (EigenResult, error) {
	n := len(original)
	for i := range original {
		for j := 0; j < i; j++ {
			if original[i][j] != original[j][i] {
				return EigenResult{}, fmt.Errorf("jacobi rotations require a symmetric matrix")
			}
		}
	}
	a, v := Clone(original), Identity(n)
	result := EigenResult{}
	for iteration := 0; iteration <= limit; iteration++ {
		p, q, largest, off := 0, 0, 0.0, 0.0
		values := make([]float64, n)
		for i := 0; i < n; i++ {
			values[i] = a[i][i]
			for j := i + 1; j < n; j++ {
				off = math.Hypot(off, a[i][j])
				if math.Abs(a[i][j]) > largest {
					p, q, largest = i, j, math.Abs(a[i][j])
				}
			}
		}
		residual := MaxDifference(Multiply(original, v), Multiply(v, Diagonal(values)))
		result.History = append(result.History, Step{iteration, values, off, residual})
		if off <= epsilon {
			result.Values, result.Vectors = values, v
			return result, nil
		}
		if iteration == limit {
			break
		}
		angle := 0.5 * math.Atan2(2*a[p][q], a[p][p]-a[q][q])
		s, c := math.Sincos(angle)
		app, aqq, apq := a[p][p], a[q][q], a[p][q]
		for i := 0; i < n; i++ {
			if i != p && i != q {
				ip, iq := a[i][p], a[i][q]
				a[i][p], a[p][i] = c*ip+s*iq, c*ip+s*iq
				a[i][q], a[q][i] = -s*ip+c*iq, -s*ip+c*iq
			}
			vp, vq := v[i][p], v[i][q]
			v[i][p], v[i][q] = c*vp+s*vq, -s*vp+c*vq
		}
		a[p][p] = c*c*app + 2*c*s*apq + s*s*aqq
		a[q][q] = s*s*app - 2*c*s*apq + c*c*aqq
		a[p][q], a[q][p] = 0, 0
	}
	return result, fmt.Errorf("jacobi rotations did not converge within %d steps", limit)
}
