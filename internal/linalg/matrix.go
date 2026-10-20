package linalg

import "math"

func Zeros(rows, cols int) Matrix {
	a := make(Matrix, rows)
	for i := range a {
		a[i] = make([]float64, cols)
	}
	return a
}

func Clone(a Matrix) Matrix {
	b := make(Matrix, len(a))
	for i := range a {
		b[i] = append([]float64(nil), a[i]...)
	}
	return b
}

func Transpose(a Matrix) Matrix {
	b := Zeros(len(a[0]), len(a))
	for i := range a {
		for j := range a[i] {
			b[j][i] = a[i][j]
		}
	}
	return b
}

func Diagonal(x []float64) Matrix {
	a := Zeros(len(x), len(x))
	for i := range x {
		a[i][i] = x[i]
	}
	return a
}

func Norm(x []float64) float64 {
	v := 0.0
	for _, x := range x {
		v = math.Max(v, math.Abs(x))
	}
	return v
}

func Distance(a, b []float64) float64 {
	v := 0.0
	for i := range a {
		v = math.Max(v, math.Abs(a[i]-b[i]))
	}
	return v
}

func NormInf(a Matrix) float64 {
	v := 0.0
	for _, row := range a {
		sum := 0.0
		for _, x := range row {
			sum += math.Abs(x)
		}
		v = math.Max(v, sum)
	}
	return v
}

type Step struct {
	Iteration int
	X         []float64
	Estimate  float64
	Residual  float64
}
