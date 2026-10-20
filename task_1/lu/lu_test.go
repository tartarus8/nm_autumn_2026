package main

import (
	"math"
	"nm/internal/linalg"
	"testing"
)

func TestLU(t *testing.T) {
	tasks, err := readTasks("examples/small.json")
	if err != nil {
		t.Fatal(err)
	}
	a, b := tasks[0].Matrix, tasks[0].RHS
	f, err := linalg.FactorLU(a, false)
	if err != nil {
		t.Fatal(err)
	}
	x := linalg.SolveLU(f, b, false)
	if !(linalg.Distance(x, []float64{1, 2, 3}) < 1e-8) {
		t.Fatalf("x = %v, want [1 2 3]", x)
	}
	if !(math.Abs(linalg.DeterminantLU(f)+3) < 1e-8) {
		t.Fatal("determinant should be -3")
	}
	if !(linalg.MaxDifference(linalg.Multiply(f.L, f.U), linalg.Multiply(f.P, a)) < 1e-8) {
		t.Fatal("LU differs from PA")
	}
	inverse := linalg.InverseLU(f, false)
	if !(linalg.MaxDifference(linalg.Multiply(a, inverse), linalg.Identity(len(a))) < 1e-8) {
		t.Fatal("A * inverse differs from I")
	}
}
