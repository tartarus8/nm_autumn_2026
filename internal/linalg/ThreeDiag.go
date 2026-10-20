package linalg

import "fmt"

func ThreeDiag(lower, diagonal, upper, rhs []float64) (x, p, q []float64, err error) {
	n := len(diagonal)
	p, q, x = make([]float64, n), make([]float64, n), make([]float64, n)
	for i := 0; i < n; i++ {
		denominator, value := diagonal[i], rhs[i]
		if i > 0 {
			denominator += lower[i-1] * p[i-1]
			value -= lower[i-1] * q[i-1]
		}
		if denominator == 0 {
			return nil, nil, nil, fmt.Errorf("ThreeDiag: zero denominator at row %d", i+1)
		}
		if i < n-1 {
			p[i] = -upper[i] / denominator
		}
		q[i] = value / denominator
	}
	x[n-1] = q[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = p[i]*x[i+1] + q[i]
	}
	return x, p, q, nil
}

func ThreeDiagProduct(lower, diagonal, upper, x []float64) []float64 {
	b := make([]float64, len(x))
	for i := range x {
		b[i] = diagonal[i] * x[i]
		if i > 0 {
			b[i] += lower[i-1] * x[i-1]
		}
		if i+1 < len(x) {
			b[i] += upper[i] * x[i+1]
		}
	}
	return b
}
