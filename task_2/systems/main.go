package main

import (
	"fmt"
	"nm/internal/labio"
	"nm/internal/linalg"
)

func run(extraInfo bool) error {
	t := task
	if err := t.Check(); err != nil {
		return err
	}
	fmt.Println("Lab 2.2, provisional variant 1: (x1^2 + 4)*x2 = 8; (x1 - 1)^2 + (x2 - 1)^2 = 4")
	fmt.Printf("epsilon = %g; initial = %v; reference = %.15g\n", t.Epsilon, t.Initial, t.Reference)
	if extraInfo {
		fmt.Printf("The curves in initial_guess.png intersect near (2.96, 0.63); choose x0 = %v.\n", t.Initial)
		fmt.Println("phi(x) = [1 + sqrt(4 - (x2-1)^2), 8/(x1^2+4)]")
		fmt.Printf("phi maps the box %v .. %v into itself; ||phi'(x)||inf <= q = %g.\n", t.Lower, t.Upper, t.Q)
		fmt.Println("Newton: solve J(x)*delta = F(x) using LU, then x_next = x - delta.")
	}
	counts := []int{}
	for _, method := range []string{"iteration", "newton"} {
		var history []linalg.Step
		var err error
		metric := "error bound"
		if method == "iteration" {
			history, err = iteration(t)
		} else {
			history, err = newton(t)
			metric = "step size"
		}
		if err != nil {
			return err
		}
		last := history[len(history)-1]
		for _, value := range last.X {
			if value <= 0 {
				return fmt.Errorf("%s found a non-positive solution", method)
			}
		}
		counts = append(counts, last.Iteration)
		fmt.Printf("%s: x = %.15g; iterations = %d\n", method, last.X, last.Iteration)
		fmt.Printf("%s = %.3e; residual = %.3e; reference error = %.3e\n", metric, last.Estimate, last.Residual, linalg.Distance(last.X, t.Reference))
		if extraInfo {
			labio.PrintHistory(method+" (estimate: "+metric+")", history, t.Reference)
		}
	}
	fmt.Printf("Simple iteration: %d iterations; Newton: %d iterations.\n", counts[0], counts[1])
	return nil
}

func main() { labio.Main(run) }
