package nonlinear

import (
	"math"
	"nm/internal/labio"
	"nm/internal/linalg"
	"testing"
)

func TestExpressions(t *testing.T) {
	e, err := Parse("pow(x1,2) + sqrt(x2) - sin(pi/2) + exp(0)", 2)
	if err != nil {
		t.Fatal(err)
	}
	if value := e([]float64{3, 4}); math.Abs(value-11) > 1e-12 {
		t.Fatalf("got %g, want 11", value)
	}
	for _, source := range []string{"x3+1", "x^2", "os.Exit(1)", "sqrt(1,2)"} {
		if _, err := Parse(source, 2); err == nil {
			t.Fatalf("accepted unsupported expression %s", source)
		}
	}
}

func TestConfiguredProblemsAndDerivatives(t *testing.T) {
	for _, file := range []string{"../../task_2/scalar/tests.json", "../../task_2/systems/tests.json"} {
		configs, err := labio.ReadTasks[Config](file)
		if err != nil {
			t.Fatal(err)
		}
		for _, c := range configs {
			t.Run(c.ID, func(t *testing.T) {
				p, err := c.Build()
				if err != nil {
					t.Fatal(err)
				}
				if linalg.Norm(p.F(c.Reference)) > 1e-12 || linalg.Distance(p.Phi(c.Reference), c.Reference) > 1e-12 {
					t.Fatal("reference is not a root and fixed point")
				}
				for _, pair := range []struct {
					f Vector
					j Jacobian
				}{{p.F, p.Jacobian}, {p.Phi, p.PhiJacobian}} {
					jac := pair.j(c.Initial)
					for j := range c.Initial {
						left, right := append([]float64(nil), c.Initial...), append([]float64(nil), c.Initial...)
						left[j] -= 1e-6
						right[j] += 1e-6
						fl, fr := pair.f(left), pair.f(right)
						for i := range jac {
							if math.Abs(jac[i][j]-(fr[i]-fl[i])/2e-6) > 1e-7 {
								t.Fatalf("incorrect derivative at (%d,%d)", i, j)
							}
						}
					}
				}
				for _, newton := range []bool{false, true} {
					history, err := Solve(p, c.Epsilon, c.MaxIterations, newton)
					if err != nil {
						t.Fatal(err)
					}
					last := history[len(history)-1]
					if linalg.Distance(last.X, c.Reference) > c.Epsilon {
						t.Fatalf("solution outside requested accuracy: %v", last.X)
					}
					if !newton {
						for _, step := range history {
							if linalg.Distance(step.X, c.Reference) > step.Estimate+1e-14 {
								t.Fatal("incorrect contraction error bound")
							}
						}
					}
				}
			})
		}
	}
}

func TestNonlinearFailures(t *testing.T) {
	configs, err := labio.ReadTasks[Config]("../../task_2/scalar/tests.json")
	if err != nil {
		t.Fatal(err)
	}
	p, err := configs[0].Build()
	if err != nil {
		t.Fatal(err)
	}
	if _, err := Solve(p, 1e-12, 1, false); err == nil {
		t.Fatal("expected iteration-limit error")
	}
	p.Q = 0.01
	if _, err := Solve(p, 1e-12, 100, false); err == nil {
		t.Fatal("expected invalid bound error")
	}
	p.Initial = []float64{-2}
	p.Lower = []float64{-3}
	if _, err := Solve(p, 1e-12, 100, false); err == nil {
		t.Fatal("expected domain error")
	}
}
