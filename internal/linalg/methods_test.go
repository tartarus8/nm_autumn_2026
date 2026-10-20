package linalg

import (
	"sort"
	"testing"
)

func check(t *testing.T, got, want Matrix) {
	t.Helper()
	if !(MaxDifference(got, want) < 1e-8) {
		t.Fatalf("got %v, want %v", got, want)
	}
}

func TestThreeDiag(t *testing.T) {
	// 2*x1 - x2 = 0; -x1 + 2*x2 - x3 = 0; -x2 + 2*x3 = 4.
	x, _, _, err := ThreeDiag([]float64{-1, -1}, []float64{2, 2, 2}, []float64{-1, -1}, []float64{0, 0, 4})
	if err != nil {
		t.Fatal(err)
	}
	check(t, Matrix{x}, Matrix{{1, 2, 3}})
}

func TestIterations(t *testing.T) {
	a, b := Matrix{{4, 1}, {1, 3}}, []float64{6, 7}
	for _, seidel := range []bool{false, true} {
		r, err := LinearIterations(a, b, 1e-10, 100, seidel)
		if err != nil {
			t.Fatal(err)
		}
		check(t, Matrix{r.X}, Matrix{{1, 2}})
	}
}

func TestRotations(t *testing.T) {
	a := Matrix{{2, 1}, {1, 2}}
	r, err := Rotations(a, 1e-10, 100)
	if err != nil {
		t.Fatal(err)
	}
	check(t, Multiply(a, r.Vectors), Multiply(r.Vectors, Diagonal(r.Values)))
	check(t, Multiply(Transpose(r.Vectors), r.Vectors), Identity(2))
	sort.Float64s(r.Values)
	check(t, Matrix{r.Values}, Matrix{{1, 3}})
}

func TestQR(t *testing.T) {
	a := Matrix{{0, -1}, {1, 0}}
	q, r := HouseholderQR(a)
	check(t, Multiply(q, r), a)
	check(t, Multiply(Transpose(q), q), Identity(2))
	result, err := QREigenvalues(a, 1e-10, 100)
	if err != nil {
		t.Fatal(err)
	}
	// The eigenvalues have real part 0 and imaginary parts -1 and +1.
	v := result.Values
	check(t, Matrix{{real(v[0]), real(v[1])}}, Matrix{{0, 0}})
	parts := []float64{imag(v[0]), imag(v[1])}
	sort.Float64s(parts)
	check(t, Matrix{parts}, Matrix{{-1, 1}})
}
