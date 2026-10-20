package linalg

import (
	"fmt"
	"math"
	"math/cmplx"
)

func reflector(x []float64) []float64 {
	norm := 0.0
	for _, value := range x {
		norm = math.Hypot(norm, value)
	}
	if norm == 0 {
		return nil
	}
	v := append([]float64(nil), x...)
	for i := range v {
		v[i] /= norm
	}
	v[0] += math.Copysign(1, x[0])
	norm = 0
	for _, value := range v {
		norm = math.Hypot(norm, value)
	}
	for i := range v {
		v[i] /= norm
	}
	return v
}

func reflectLeft(a Matrix, v []float64, row, col int) {
	for j := col; j < len(a[0]); j++ {
		dot := 0.0
		for i := range v {
			dot += v[i] * a[row+i][j]
		}
		for i := range v {
			a[row+i][j] -= 2 * v[i] * dot
		}
	}
}

func reflectRight(a Matrix, v []float64, col int) {
	for i := range a {
		dot := 0.0
		for j := range v {
			dot += a[i][col+j] * v[j]
		}
		for j := range v {
			a[i][col+j] -= 2 * dot * v[j]
		}
	}
}

func HouseholderQR(a Matrix) (q, r Matrix) {
	n := len(a)
	q, r = Identity(n), Clone(a)
	for k := 0; k < n-1; k++ {
		x := make([]float64, n-k)
		for i := range x {
			x[i] = r[k+i][k]
		}
		v := reflector(x)
		if v == nil {
			continue
		}
		reflectLeft(r, v, k, k)
		reflectRight(q, v, k)
		for i := k + 1; i < n; i++ {
			r[i][k] = 0
		}
	}
	return
}

func hessenberg(a Matrix) Matrix {
	h := Clone(a)
	for k := 0; k < len(a)-2; k++ {
		x := make([]float64, len(a)-k-1)
		for i := range x {
			x[i] = h[k+1+i][k]
		}
		v := reflector(x)
		if v == nil {
			continue
		}
		reflectLeft(h, v, k+1, k)
		reflectRight(h, v, k+1)
		for i := k + 2; i < len(a); i++ {
			h[i][k] = 0
		}
	}
	return h
}

func eigen2(a, b, c, d float64) (complex128, complex128) {
	scale := math.Max(math.Max(math.Abs(a), math.Abs(b)), math.Max(math.Abs(c), math.Abs(d)))
	if scale == 0 {
		return 0, 0
	}
	a, b, c, d = a/scale, b/scale, c/scale, d/scale
	center := complex((a+d)/2, 0)
	root := cmplx.Sqrt(complex((a-d)*(a-d)/4+b*c, 0))
	return (center + root) * complex(scale, 0), (center - root) * complex(scale, 0)
}

type QRStep struct {
	Iteration, Size int
	Shift           string
	Subdiagonal     float64
}
type QRResult struct {
	Values  []complex128
	History []QRStep
}

func QREigenvalues(a Matrix, epsilon float64, limit int) (QRResult, error) {
	h := hessenberg(a)
	result := QRResult{}
	m, steps, stalled := len(a), 0, 0
	for m > 0 {
		low := m - 1
		for low > 0 {
			scale := math.Abs(h[low-1][low-1]) + math.Abs(h[low][low])
			if math.Abs(h[low][low-1]) <= epsilon*math.Max(1, scale) {
				h[low][low-1] = 0
				break
			}
			low--
		}
		size := m - low
		if size == 1 {
			result.Values = append(result.Values, complex(h[m-1][m-1], 0))
			m--
			stalled = 0
			continue
		}
		if size == 2 {
			first, second := eigen2(h[low][low], h[low][low+1], h[low+1][low], h[low+1][low+1])
			result.Values = append(result.Values, first, second)
			m -= 2
			stalled = 0
			continue
		}
		if steps >= limit {
			return result, fmt.Errorf("QR did not converge within %d iterations", limit)
		}
		block := Zeros(size, size)
		for i := range block {
			copy(block[i], h[low+i][low:m])
		}
		scale := NormInf(block)
		for i := range block {
			for j := range block {
				block[i][j] /= scale
			}
		}
		last := size - 1
		first, second := eigen2(block[last-1][last-1], block[last-1][last], block[last][last-1], block[last][last])
		var next Matrix
		kind := "single"
		if imag(first) == 0 || stalled > 0 && stalled%30 == 0 {
			mu := real(first)
			if math.Abs(real(second)-block[last][last]) < math.Abs(mu-block[last][last]) {
				mu = real(second)
			}
			if stalled > 0 && stalled%30 == 0 {
				mu = block[last][last] + math.Abs(block[last][last-1])
				kind = "exceptional"
			}
			shifted := Clone(block)
			for i := range shifted {
				shifted[i][i] -= mu
			}
			q, r := HouseholderQR(shifted)
			next = Multiply(r, q)
			for i := range next {
				next[i][i] += mu
			}
		} else {
			kind = "double"
			trace := block[last-1][last-1] + block[last][last]
			det := block[last-1][last-1]*block[last][last] - block[last-1][last]*block[last][last-1]
			polynomial := Multiply(block, block)
			for i := range polynomial {
				for j := range polynomial {
					polynomial[i][j] -= trace * block[i][j]
				}
				polynomial[i][i] += det
			}
			q, _ := HouseholderQR(polynomial)
			next = Multiply(Transpose(q), Multiply(block, q))
		}
		for i := range next {
			for j := range next {
				if i > j+1 {
					next[i][j] = 0
				}
				h[low+i][low+j] = next[i][j] * scale
			}
		}
		steps++
		stalled++
		measure := math.Abs(h[m-1][m-2])
		if math.IsNaN(measure) || math.IsInf(measure, 0) {
			return result, fmt.Errorf("non-finite QR iteration")
		}
		result.History = append(result.History, QRStep{steps, size, kind, measure})
	}
	return result, nil
}
