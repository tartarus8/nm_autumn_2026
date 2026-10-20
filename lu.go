package main

import (
	"fmt"
	"math"
)

type Matrix [][]float64

type Factors struct {
	L, U, P Matrix
	Swaps   int
}

func identity(n int) Matrix {
	a := make(Matrix, n)
	for i := range a {
		a[i] = make([]float64, n)
		a[i][i] = 1
	}
	return a
}

// factorLU uses partial pivoting to construct P*A = L*U.
// The input is a nonempty square matrix and is left unchanged.
func factorLU(a Matrix, extraInfo bool) (Factors, error) {
	n := len(a)
	f := Factors{L: identity(n), U: make(Matrix, n), P: identity(n)}
	for i := range a {
		f.U[i] = append([]float64(nil), a[i]...)
	}

	for k := 0; k < n; k++ {
		pivot := k
		for i := k + 1; i < n; i++ {
			if math.Abs(f.U[i][k]) > math.Abs(f.U[pivot][k]) {
				pivot = i
			}
		}
		if f.U[pivot][k] == 0 {
			return Factors{}, fmt.Errorf("singular matrix: zero pivot in column %d", k+1)
		}
		if extraInfo {
			fmt.Printf("Column %d: choose row %d, pivot = %.8g\n", k+1, pivot+1, f.U[pivot][k])
		}

		if pivot != k {
			f.U[k], f.U[pivot] = f.U[pivot], f.U[k]
			f.P[k], f.P[pivot] = f.P[pivot], f.P[k]
			// Only completed columns of L contain elimination history.
			for j := 0; j < k; j++ {
				f.L[k][j], f.L[pivot][j] = f.L[pivot][j], f.L[k][j]
			}
			f.Swaps++
			if extraInfo {
				fmt.Printf("Swap rows %d and %d in U and P; swap the first %d entries of these rows in L\n",
					k+1, pivot+1, k)
			}
		}

		for i := k + 1; i < n; i++ {
			f.L[i][k] = f.U[i][k] / f.U[k][k]
			if extraInfo {
				fmt.Printf("L[%d,%d] = %.8g / %.8g = %.8g; R%d <- R%d - (%.8g)*R%d\n",
					i+1, k+1, f.U[i][k], f.U[k][k], f.L[i][k], i+1, i+1, f.L[i][k], k+1)
			}
			f.U[i][k] = 0
			for j := k + 1; j < n; j++ {
				f.U[i][j] -= f.L[i][k] * f.U[k][j]
			}
		}
		if extraInfo {
			printMatrix("L after this column", f.L)
			printMatrix("U after this column", f.U)
			printMatrix("P after this column", f.P)
		}
	}
	return f, nil
}

func matVec(a Matrix, x []float64) []float64 {
	result := make([]float64, len(a))
	for i := range a {
		for j := range a[i] {
			result[i] += a[i][j] * x[j]
		}
	}
	return result
}

func solveLU(f Factors, b []float64, extraInfo bool) []float64 {
	n := len(b)
	pb := matVec(f.P, b)
	if extraInfo {
		fmt.Printf("Pb = %.10g\nForward substitution: Ly = Pb\n", pb)
	}

	y := make([]float64, n)
	for i := 0; i < n; i++ {
		y[i] = pb[i]
		for j := 0; j < i; j++ {
			y[i] -= f.L[i][j] * y[j]
		}
		if extraInfo {
			fmt.Printf("y[%d] = %.8g\n", i+1, y[i])
		}
	}

	if extraInfo {
		fmt.Println("Backward substitution: Ux = y")
	}
	x := make([]float64, n)
	for i := n - 1; i >= 0; i-- {
		x[i] = y[i]
		for j := i + 1; j < n; j++ {
			x[i] -= f.U[i][j] * x[j]
		}
		if extraInfo {
			fmt.Printf("x[%d] = %.10g / %.10g", i+1, x[i], f.U[i][i])
		}
		x[i] /= f.U[i][i]
		if extraInfo {
			fmt.Printf(" = %.10g\n", x[i])
		}
	}
	return x
}

func determinantLU(f Factors) float64 {
	det := 1.0
	if f.Swaps%2 != 0 {
		det = -1
	}
	for i := range f.U {
		det *= f.U[i][i]
	}
	return det
}

func inverseLU(f Factors, extraInfo bool) Matrix {
	n := len(f.U)
	inverse := make(Matrix, n)
	for i := range inverse {
		inverse[i] = make([]float64, n)
	}
	for j := 0; j < n; j++ {
		e := make([]float64, n)
		e[j] = 1
		if extraInfo {
			fmt.Printf("\nInverse column %d: solve A*z = e%d, e%d = %v using the same factors\n",
				j+1, j+1, j+1, e)
		}
		column := solveLU(f, e, extraInfo)
		for i := 0; i < n; i++ {
			inverse[i][j] = column[i]
		}
	}
	return inverse
}

func multiply(a, b Matrix) Matrix {
	c := make(Matrix, len(a))
	for i := range a {
		c[i] = make([]float64, len(b[0]))
		for k := range b {
			for j := range b[k] {
				c[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return c
}

func maxDifference(a, b Matrix) float64 {
	difference := 0.0
	for i := range a {
		for j := range a[i] {
			difference = math.Max(difference, math.Abs(a[i][j]-b[i][j]))
		}
	}
	return difference
}
